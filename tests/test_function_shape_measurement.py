"""Hand-calculated inventory independent of complexity tools."""

from tools.measure_function_shape import function_shapes


def test_function_spans_and_arguments_do_not_double_count_classes():
    source = """class Example:
    def method(self, a, /, b=1, *args, named=2, **kwargs):
        def nested(value):
            return value
        return nested(a)

async def operation():
    pass
"""
    rows = function_shapes(source, "sample.py")
    assert [(r["function"], r["lines"], r["parameters"]) for r in rows] == [
        ("Example.method", 4, 6),
        ("Example.method.nested", 2, 1),
        ("operation", 2, 0),
    ]
    assert [(r["line"], r["end_line"]) for r in rows] == [(2, 5), (3, 4), (7, 8)]
