# EXPERIMENTAL
import operator as operators

def value(expr, **kwargs):
    if type(expr) in (Var, Expr):
        return expr.eval(**kwargs)

    return expr

class Expr:
    def __init__(self, operator, operands):
        self.operator = operator
        self.operands = operands

    def __pos__(self):
        return self

    def __neg__(self):
        return -1 * self

    def __add__(self, other):
        return Expr("+", [self, other])

    def __radd__(self, other):
        return Expr("+", [other, self])

    def __sub__(self, other):
        return Expr("-", [self, other])

    def __rsub__(self, other):
        return Expr("-", [other, self])

    def __mul__(self, other):
        return Expr("*", [self, other])

    def __rmul__(self, other):
        return Expr("*", [other, self])

    def __truediv__(self, other):
        return Expr("/", [self, other])

    def __rtruediv__(self, other):
        return Expr("/", [other, self])

    def __pow__(self, other):
        return Expr("^", [self, other])

    def __rpow__(self, other):
        return Expr("^", [other, self])

    def eval(self, **kwargs):
        operator_map = {
            "+": operators.add,
            "-": operators.sub,
            "*": operators.mul,
            "/": operators.truediv,
            "^": operators.pow
        }

        operator = operator_map[self.operator]
        acc = value(self.operands[0], **kwargs)
        for other in self.operands[1:]:
            acc = operator(acc, value(other, **kwargs))

        return acc

    def __repr__(self):
        if len(self.operands) == 1:
            return "(" + self.operator + repr(self.operands[0]) + ")"

        return "(" + self.operator.join(map(repr, self.operands)) + ")"

class Var(Expr):
    def __init__(self, name):
        self.name = name

    def eval(self, **kwargs):
        return kwargs[self.name]

    def __repr__(self):
        return self.name
