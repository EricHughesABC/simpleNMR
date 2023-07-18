from collections.abc import Iterable

def is_iterable(obj):
    return isinstance(obj, Iterable)


def stringify_vals(vals, separator=", ") -> str:
    if not is_iterable(vals):
        vals_set = {vals}
    elif isinstance(vals, str):
        vals_set = {x.strip() for x in vals.split(",")}
    else:
        vals_set = set(vals)
        # rejoin the set into a string
    return separator.join([f"{str(v)}" for v in vals_set])