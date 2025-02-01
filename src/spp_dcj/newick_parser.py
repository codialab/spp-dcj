from collections import deque
from io import BytesIO, StringIO
import re

IS_FLOAT = re.compile(r'^\s*(?:-?((\d*\.)?\d+)|((\d*\.\d+|\d+)(e|E)(\+|-)?\d+))\s*$')
token_pat = re.compile(r'(?:\(|\)|,|\:|;|[^():;,]+)')

class ParseType:

    TREE='TREE'
    BRANCH='BRANCH'
    LEAF='LEAF'
    STRING='STRING'
    FLOAT='FLOAT'

def parse_tree_iterator(data):
    if type(data) == str:
        data = StringIO(data)
    cursor = data.read(1)
    program = StringIO()
    while cursor:
        program.write(cursor)
        next_c = data.read(1)
        if cursor==';' or not next_c:
            t = __parse_tree__(program.getvalue())
            if t:
                yield t[0]
            program = StringIO()
        cursor = next_c

def parse_tree(data):
    res = []
    if type(data) == str:
        data = StringIO(data)
    cursor = data.read(1)
    program = StringIO()
    while cursor:
        program.write(cursor)
        next_c = data.read(1)
        if cursor==';' or not next_c:
            res.extend(__parse_tree__(program.getvalue()))
            program = StringIO()
        cursor = next_c
    if len(res) == 1:
        return res[0]
    return res

def __parse_tree__(program):
    res = []
    stack = [ParseType.TREE]
    for token in token_pat.findall(program):
        prev = stack[-1]
        if token == '(':
            if prev != ParseType.LEAF and prev != ParseType.TREE:
                raise SyntaxError('could not parse tree')
            # remove faulty prediction
            if prev == ParseType.LEAF:
                stack.pop()
            stack.append(ParseType.BRANCH)
            stack.append(ParseType.LEAF)
        elif token == ')':
            prev = stack.pop()
            subtrees = deque()
            while(prev != ParseType.BRANCH):
                subtrees.appendleft(prev)
                if not len(stack):
                    raise SyntaxError('too many closing brackets..')
                prev = stack.pop()
            stack.append(Branch(subtrees))
        elif token == ',':
            # predict leaf
            stack.append(ParseType.LEAF)
        elif token == ':':
            stack.append(ParseType.FLOAT)
        elif token == ';':
            if len(stack) > 1:
                prev = stack.pop()
                subtrees = deque()
                while(prev != ParseType.TREE):
                    if type(prev) != Leaf and type(prev) != Branch:
                        raise SyntaxError('unexpected end of tree')
                        continue
                    subtrees.appendleft(prev)
                    prev = stack.pop()
                if len(subtrees) > 1:
                    subtrees = [Branch(subtrees)]
                stack.append(Tree(subtrees[0]))
            res.append(stack.pop())
            if len(stack) != 0:
                raise SyntaxError('could not reduce tree...')
            stack.append(ParseType.TREE)
        elif not token.strip():
            # ignore whitespaces
            continue
        elif prev == ParseType.LEAF:
            # remove LEAF flag
            stack.pop()
            stack.append(Leaf(token.strip()))
        elif prev == ParseType.FLOAT:
            if IS_FLOAT.match(token):
                # remove FLOAT flag
                stack.pop()
                stack[-1].length = float(token)
            else:
                raise SyntaxError('%s was expected to be a number' %token)
        elif type(prev) == Branch:
            prev.label = token.strip()
        elif prev == ParseType.TREE:
            stack.append(Leaf(token.strip()))
        else:
            raise SyntaxError('unknown operator')
    if len(stack) > 1 or (stack and stack[0] != ParseType.TREE):
        raise SyntaxError('Error while parsing tree')
    return res

class Leaf(object):
    
    def __init__(self, label, length=None):
        self.label = label
        self.length = length
        self.ancestor = None

    def setAncestor(self, anc):
        self.ancestor = anc

    def getAncestor(self):
        return self.ancestor

    def __str__(self):
        return '%s%s' % (self.label, self.length != None and ':%.5f' %self.length or '')

class Branch(object):

    def __init__(self, subtrees, label=None, length=None):
        self.subtrees = list()
        self.ancestor = None
        for t in subtrees:
            if type(t) != Branch and type(t) != Leaf:
                raise TypeError('%s is not of type Leaf or Branch but of %s' %(t, type(t)))
            t.setAncestor(self)
            self.subtrees.append(t)
        self.label = label or ''
        self.length = length

    def setAncestor(self, anc):
        self.ancestor = anc

    def getAncestor(self):
        return self.ancestor

    def __str__(self):
        return '(%s)%s%s' % (','.join(map(str, self.subtrees)),self.label != None
                and self.label or '', self.length != None and ':%.5f' %self.length or '')

    def getLeaves(self):
        res = []
        for t in self.subtrees:
            if type(t) == Leaf:
                res.append(t)
            elif type(t) == Branch:
                res.extend(t.getLeaves())
        return res

    def getNodes(self):
        res = []
        for t in self.subtrees:
            res.append(t)
            if type(t) == Branch:
                res.extend(t.getNodes())
        return res


class Tree(object):
    def __init__(self, subtree):
        self.subtree = subtree
    def __str__(self):
        return str(self.subtree) + ';' 
    def getLeaves(self):
        return self.subtree.getLeaves()
    def getNodes(self):
        return self.subtree.getNodes()

