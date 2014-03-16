
# Metalin Copyright (C) 2014 2-Complex
# 
# This software is provided 'as-is', without any express or implied
# warranty.  In no event will the authors be held liable for any damages
# arising from the use of this software.
# 
# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:
# 
# 1. The origin of this software must not be misrepresented; you must not
#  claim that you wrote the original software. If you use this software
#  in a product, an acknowledgment in the product documentation would be
#  appreciated but is not required.
# 2. Altered source versions must be plainly marked as such, and must not be
#  misrepresented as being the original software.
# 3. This notice may not be removed or altered from any source distribution.


NORMALIZE_FUNCTION = False
MAXIMUM_DIMENSION = 5
SCALAR_TYPE = 'double'
VARIABLE_NAMES = 'xyzwt'


def do(t, j, vars):
    l = []
    i = 0
    for v in vars:
        l += [t.replace('x', v).replace('_i', str(i))]
        i+=1
    return j.join(l)

def vector_space_ops(scalar, classname, vars):
    '''Generates code for vector space operators addition subtraction, scalar multiplication.'''
    s = ''
    s+= '\tinline ' + classname + '& operator= (const ' + classname + '& v) {' + do('x=v.x; ', '', vars) + 'return *this;}\n'
    s+= '\tinline bool operator== (const ' + classname + '& v) const {return ' + do('x==v.x', ' && ', vars) + ';}\n'
    s+= '\tinline bool operator!= (const ' + classname + '& v) const {return ' + do('x!=v.x', ' || ', vars) + ';}\n'
    s+= '\tinline void set(' + do('' + scalar + ' in_x', ', ', vars) + ') {' + do('x=in_x;' ,' ', vars) + '}\n'
    s+= '\tinline void set(const ' + scalar + '* v) {' + do('x=v[_i];' ,' ', vars) + '}\n'
    s+= '\tinline ' + classname + ' operator- () const {return ' + classname + '(' + do('-x', ', ', vars) + ');}\n'
    s+= '\tinline ' + classname + ' operator+ (const ' + classname + '& v) const {return ' + classname + '(' + do('x+v.x', ', ', vars) + ');}\n'
    s+= '\tinline ' + classname + ' operator- (const ' + classname + '& v) const {return ' + classname + '(' + do('x-v.x', ', ', vars) + ');}\n'
    s+= '\tinline ' + classname + ' operator* (' + scalar + ' k) const {return ' + classname + '(' + do('x*k', ', ', vars) + ');}\n'
    s+= '\tinline ' + classname + ' operator/ (' + scalar + ' k) const {return ' + classname + '(' + do('x/k', ', ', vars) + ');}\n'
    s+= '\tinline friend ' + classname + ' operator* (' + scalar + ' k, const ' + classname + '& v) {return ' + classname + '(' + do('k*v.x', ', ', vars) + ');}\n'
    s+= '\tinline friend ' + classname + ' operator/ (' + scalar + ' k, const ' + classname + '& v) {return ' + classname + '(' + do('k/v.x', ', ', vars) + ');}\n'
    s+= '\tinline ' + classname + '& operator+= (const ' + classname + '& v) {' + do('x+=v.x; ', '', vars) + 'return *this;}\n'
    s+= '\tinline ' + classname + '& operator-= (const ' + classname + '& v) {' + do('x-=v.x; ', '', vars) + 'return *this;}\n'
    s+= '\tinline ' + classname + '& operator*= (' + scalar + ' k) {' + do('x*=k; ', '', vars) + 'return *this;}\n'
    s+= '\tinline ' + classname + '& operator/= (' + scalar + ' k) {' + do('x/=k; ', '', vars) + 'return *this;}\n'
    return s

def vector_code(scalar, classbase, vars):
    '''Generates code for a vector class.'''
    n = str(len(vars))
    classname = classbase + str(n)
    s = ''
    s+= 'class ' + classname + ';\n'
    s+= 'class ' + classname + ' {\n'
    s+= 'public:\n'
    s+= '\t' + scalar + ' ' + ', '.join(vars) + ';\n'
    s+= '\t' + classname + '() : ' + do('x(0.0)', ', ', vars) + ' {}\n'
    s+= constructors(scalar, '\t', classbase, vars)
    s+= '\t\n'
    s+= vector_space_ops(scalar, classname, vars)
    s+= '\t\n'
    s+= '\tinline ' + classname + ' operator* (const ' + classname + '& v) const {return ' + classname + '(' + do('x*v.x', ', ', vars) + ');}\n'
    s+= '\tinline ' + classname + ' operator/ (const ' + classname + '& v) const {return ' + classname + '(' + do('x/v.x', ', ', vars) + ');}\n'
    s+= '\tinline ' + classname + '& operator*= (const ' + classname + '& v) {' + do('x*=v.x; ', '', vars) + 'return *this;}\n'
    s+= '\tinline ' + classname + '& operator/= (const ' + classname + '& v) {' + do('x/=v.x; ', '', vars) + 'return *this;}\n'
    s+= '\t\n'
    s+= '\t' + scalar + ' dot(const ' + classname + '& v) const {return ' + do('x*v.x', ' + ', vars) + ';}\n'
    s+= '\t' + scalar + ' magSquared() const {return ' + do('x*x', ' + ', vars) + ';}\n'
    s+= '\t' + scalar + ' mag() const {return sqrt(' + do('x*x', ' + ', vars) + ');}\n'
    if n=='3':
        s+= '\t' + classname + ' cross(const ' + classname + '& v) const {return ' + classname + '(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);}\n'
    s+= '\t\n'
    s+= '\tconst ' + scalar + '* ptr() const {return &(' + vars[0] + ');}\n'
    s+= '\tvoid print() const {printf("(' + do('%f', ', ', vars) + ')", ' + do('x', ', ', vars) + ');}\n'
    s+= '\tvoid display() const {print(); printf("\\n");}\n'
    s+= '};\n'
    if NORMALIZE_FUNCTION:
        s+= classname + ' normalize(const ' + classname + '& v) {return v/v.mag();}\n'
    return s

def ways(n):
    if n==0:
        return [[]]
    if n==1:
        return [[1]]
    L = []
    for i in range(1, n+1):
        w = ways(n-i)
        for l in w:
            L += [[i] + l]
    return L

def vars_sep(vars, l):
    r = []
    j=0
    for i in l:
        r += [vars[j:(j+i)]]
        j += i
    return r

def constructors(scalar, indent, classbase, vars):
    result = ''
    n = len(vars)
    classname = classbase + str(n)
    for w in ways(n):
        l = vars_sep(vars, w)
        arg_list = []
        copy_list = []
        for s in l:
            arg_list += [[scalar, 'const ' + classbase + str(len(s)) + '&'][len(s) > 1] + ' ' + s]
            if len(s)==1:
                copy_list += [s + '(' + s + ')']
            else:
                i = 0
                for c in s:
                    copy_list += [c + '(' + s + '.' + vars[i] + ')']
                    i+=1
        result += indent + classname + '(' + ', '.join(arg_list) + ') : ' + ', '.join(copy_list) + ' {}\n'
    return result



def dirichlet(i,j):
    if i==j:
        return '1.0'
    return '0.0';

def identity(n):
    l = []
    for i in range(0, n):
        for j in range(0, n):
            l += ['m%d%d(%s)'%(i,j,dirichlet(i,j))]
    return ', '.join(l)

def parens(s, chars):
    t = 0
    a = ''
    for c in s:
        t += {'(':1, ')':-1}.get(c, 0)
        a += {0:c}.get(t, '')
    for c in chars:
        if c in a:
            return '(', ')'
    return '', ''


class expression:
    s = ''
    def __init__(self, s):
        self.s = s
    
    def __str__(self):
        return '' + self.s
    
    def __repr__(self):
        return 'expression(\'' + self.s + '\')'
    
    def __add__(self, other):
        return expression(self.s + ' + ' + other.s)
    
    def __sub__(self, other):
        p0, p1 = parens(other.s, '+-')
        return expression(self.s + ' - ' + p0 + other.s + p1)

    def __mul__(self, other):
        if type(other) == int or type(other) == float:
            other = expression(str(other))
        p0, p1, p2, p3 = parens(self.s, '+-') + parens(other.s, '+-')
        return expression(p0 + self.s + p1 + '*' + p2 + other.s + p3)
    
    def __rmul__(self, other):
        if type(other) == int or type(other) == float:
            other = expression(str(other))
        p0, p1, p2, p3 = parens(self.s, '+-') + parens(other.s, '+-')
        return expression(p0 + other.s + p1 + '*' + p2 + self.s + p3)

    def __div__(self, other):
        if type(other) == int or type(other) == float:
            other = expression(str(other))
        p0, p1, p2, p3 = parens(self.s, '+-') + parens(other.s, '+-')
        return expression(p0 + self.s + p1 + '/' + p2 + other.s + p3)

    def __rdiv__(self, other):
        if type(other) == int or type(other) == float:
            other = expression(str(other))
        p0, p1, p2, p3 = parens(self.s, '+-') + parens(other.s, '+-')
        return expression(p0 + other.s + p1 + '/' + p2 + self.s + p3)

    def __iadd__(self, other):
        self.s = (self + other).s
        return self

    def __isub__(self, other):
        self.s = (self - other).s
        return self

    def __imul__(self, other):
        self.s = (self * other).s
        return self

    def __idiv__(self, other):
        self.s = (self / other).s
        return self

    def __neg__(self):
        p0, p1 = parens(self.s, '+-')
        return expression('-' + p0 + self.s + p1)


def matrix(n, s='m'):
    m = []
    for i in range(0, n):
        r = []
        for j in range(0, n):
            r += [expression(s + str(i) + str(j))]
        m += [r]
    return m


def column_vector(vars):
    return [map(expression, vars)]


def row_vector(vars):
    def foo(x):
        return [expression(x)]
    return map(foo, vars)


def cofactor(m, i, j):
    c = eval(repr(m))
    c.pop(i);
    for k in range(0, len(c)):
        c[k].pop(j);
    return c


def det(m):
    if len(m) == 1:
        return m[0][0]
    else:
        r = m[0]
        d = 0
        for i in range(0, len(r)):
            k = r[i] * det(cofactor(m, 0, i))
            if i==0:
                d = k
            else:
                if i%2:
                    d -= k
                else:
                    d += k
    return d


def mul(a, b):
    m = len(a[0])  # Rows of a = rows of answer
    n = len(b)     # Columns of b = columns of answer
    r = len(a)     # Columns of a = length of dot product
    c = []
    for i in range(0, n): # i is a column number
        c += [[0.0]*m]
        for j in range(0, m): # j is a row number
            for k in range(0, r):
                x = a[k][j] * b[i][k]
                if c[i][j] == 0.0:
                    c[i][j] = x
                else:
                    c[i][j] += x
    return c


def part_inverse(m):
    n = len(m)
    if n==1:
        return [[1.0 / m[0][0]]]
    else:
        c = []
        d = det(m)
        for i in range(0, n):
            c += [[0.0]*n]
            for j in range(0, n):
                k = det(cofactor(m, j, i))
                if c[i][j]==0.0:
                    if (i+j)%2:
                        c[i][j] = -k
                    else:
                        c[i][j] = k
                else:
                    if (i+j)%2:
                        c[i][j] -= k
                    else:
                        c[i][j] += k
        return c


def matrix_vars(n):
    vars = []
    for i in range(0, n):
        for j in range(0, n):
            vars += ['m%d%d'%(i,j)]
    return vars


def matrix_mul(scalar, n, classname):
    s = ''
    s+= '\tinline ' + classname + ' operator* (const ' + classname + '& m) const {\n\t\treturn ' + classname + '(\n'
    product = mul(matrix(n), matrix(n, 'm.m'))
    for i in range(0, n):
        for j in range(0, n):
            s+= '\t\t\t' + product[i][j].s + ['', ','][i!=n-1 or j!=n-1] + '\n'
    s+= '\t\t);\n\t}\n'
    return s

def matrix_vector_mul(scalar, n, matrix_classname, vector_classname, vars):
    s = ''
    s+= '\t' + vector_classname + ' operator* (const ' + vector_classname + '& v) const {\n\t\treturn ' + vector_classname + '(\n'
    def ap(x):
        return 'v.' + x
    product = mul(matrix(n), column_vector(map(ap, vars)))
    for i in range(0, n):
        s+= '\t\t\t' + product[0][i].s + ['', ','][i!=n-1] + '\n'
    s+= '\t\t);\n\t}\n'
    return s

def vector_matrix_mul(scalar, n, matrix_classname, vector_classname, vars):
    s = ''
    s+= '\tfriend ' + vector_classname + ' operator* (const ' + vector_classname + '& v, const ' + matrix_classname + '& m) {\n\t\treturn ' + vector_classname + '(\n'
    def ap(x):
        return 'v.' + x
    product = mul(row_vector(map(ap, vars)), matrix(n, 'm.m'))
    for i in range(0, n):
        s+= '\t\t\t' + product[i][0].s + ['', ','][i!=n-1] + '\n'
    s+= '\t\t);\n\t}\n'
    return s

def matrix_from_vectors(scalar, n, matrix_classname, mars, vector_classname, vars):
    ms = []
    for i in range(0, n):
        ms += ['m'+str(i)]
    inits = []
    for i in range(0, n):
        j = 0
        for v in vars:
            inits += ['m' + str(i) + str(j) + '(m' + str(i) + '.' + v + ')']
            j+=1
    return matrix_classname + '(' + do('const ' + vector_classname + '& x', ', ', ms) + ') : ' + ', '.join(inits) + ' {}\n'

def determinant_code(scalar, n, classname):
    s = ''
    s += '\t' + scalar + ' det() const {\n'
    body = '\t\treturn ' + det(matrix(n)).s + ';\n'
    s += make_temps(scalar, '\t\t', body);
    s += '\t}\n'
    return s

def transpose_code(scalar, n, classname):
    s = ''
    s += '\t' + classname + ' transpose() const {\n'
    s += '\t\treturn ' + classname + '(\n'
    for i in range(0, n):
        s += '\t\t\t'
        for j in range(0, n):
            s += 'm' + str(j) + str(i) + ['', ', '][j!=n-1]
        s += ['', ','][i!=n-1] + '\n'
    s += '\t\t);\n'
    s += '\t}\n'
    return s

def inverse_code(scalar, n, classname):
    mi = part_inverse(matrix(n))
    s = ''
    s += '\t' + classname + ' inverse() const {\n'
    body = ''
    d = det(matrix(n)).s
    body += '\t\t' + scalar + ' ' + 'd = ' + d + ';\n'
    body += '\t\treturn ' + classname + '(\n'
    for i in range(0, n):
        for j in range(0, n):
            body += '\t\t\t(' + mi[i][j].s + ')' + ['', ','][i!=n-1 or j!=n-1] + '\n'
    body += '\t\t) / d;'
    body = make_temps(scalar, '\t\t', body)
    s += body;
    s += '\n\t}\n'
    return s

def matrix_code(scalar, n, matrix_classname, mars, vector_classname, vars):
    s = ''
    s+= 'class ' + matrix_classname + ' {\n'
    s+= 'public:\n'
    s+= '\t' + scalar + ' ' + ', '.join(mars) + ';\n'
    s+= '\t' + matrix_classname + '() : ' + identity(n) + ' {}\n'
    s+= '\t' + matrix_classname + '(' + do(scalar + ' x', ', ', mars) + ') : ' + do('x(x)', ', ', mars) + ' {}\n'
    s+= '\t' + matrix_from_vectors(scalar, n, matrix_classname, mars, vector_classname, vars)
    s+= '\t\n'
    s+= vector_space_ops(scalar, matrix_classname, mars)
    s+= matrix_mul(scalar, n, matrix_classname)
    s+= matrix_vector_mul(scalar, n, matrix_classname, vector_classname, vars)
    s+= vector_matrix_mul(scalar, n, matrix_classname, vector_classname, vars)
    s+= determinant_code(scalar, n, matrix_classname)
    s+= inverse_code(scalar, n, matrix_classname)
    s+= transpose_code(scalar, n, matrix_classname)
    s+= '\n'
    s+= '\tconst ' + scalar + '* ptr() const {return &(' + mars[0] + ');}\n'
    s+= '\tvoid print() const {printf("(' + do('%f', ', ', mars) + ')", ' + do('x', ', ', mars) + ');}\n'
    s+= '\tvoid display() const {print(); printf("\\n");}\n'
    s+= '};\n'
    return s


def best_rep(code):
    D = {}
    L = ['']
    for c in code:
        if c=='(':
            L += ['']
        for i in range(0, len(L)):
            L[i] += c
        if c==')':
            s = L.pop()
            if D.has_key(s):
                D[s] += 1
            else:
                D[s] = 1
    best_freq = 0
    best_length = 0
    best_string = ''
    for s in D:
        if D[s] > 1 and (len(s) > best_length or (len(s) == best_length and D[s] > best_freq)):
            best_freq = D[s]
            best_length = len(s)
            best_string = s
    if best_freq > 1:
        return best_string
    else:
        return ''


def make_temps(scalar, indent, s):
    temp_count = 0
    expression = best_rep(s)
    s = '/****/\n' + s
    while expression != '' and expression[:4] != 'temp':
        temp_name = 'temp_' + str(temp_count)
        where = min(s.index('/****/'), s.index(expression))
        while where > 0 and s[where-1] != '\n':
            where-=1
        s = s.replace(expression, temp_name)
        s = s[:where] + indent + scalar + ' ' + temp_name + ' = ' + expression[1:(len(expression)-1)] + ';\n' + s[where:]
        expression = best_rep(s)
        temp_count += 1
    s = s.replace('/****/\n', '')
    return s



cpplicense = """
/*
  Metalin Copyright (C) 2012 2-Complex

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/
"""

lintop = '''
#ifndef _LIN_
#define _LIN_

#include <stdio.h>
#include <math.h>
'''

linbottom = '''
#endif
'''

def lingen():
    max_dim = MAXIMUM_DIMENSION
    scalar = SCALAR_TYPE
    vars = VARIABLE_NAMES
    
    print cpplicense

    print lintop
    
    for n in range(2, max_dim+1):
        print vector_code(scalar, 'Vec', vars[:n])
    
    for n in range(2, max_dim+1):
        print matrix_code(scalar, n, 'Mat' + str(n), matrix_vars(n), 'Vec' + str(n), vars[:n])
    
    print linbottom


import sys

if __name__ == "__main__":
        lingen()


