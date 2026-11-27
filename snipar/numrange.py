import argparse
import re

def parseNumRange(string):
    """reads either a int or a range"""
    match = re.fullmatch(r"\s*(\d+)(?:\s*-\s*(\d+))?\s*", string)
    if match is None:
        raise argparse.ArgumentTypeError(f"invalid integer or range: {string}")
    start = match[1]
    if match[2] is None:
        return start
    end = int(match[2])
    if end < int(start):
        raise argparse.ArgumentTypeError(f"range must be in ascending order: {string}")
    return [str(n) for n in range(int(start), end + 1)]


class NumRangeAction(argparse.Action):
    """flattens and sorts the resulting num range. also removes duplicates"""
    def __call__(self, parser, args, values, option_string=None):
        # Flatten the input values: extend the result with list elements or append string elements.
        result = []
        for v in values:
            if isinstance(v, list):
                result.extend(v)
            elif isinstance(v, str):
                result.append(v)
        
        # Convert the string values to integers, remove duplicates with a set, and sort them.
        unique_ints = sorted({int(x) for x in result})
        
        # Convert the integers back to strings.
        final_result = [str(x) for x in unique_ints]
        
        setattr(args, self.dest, final_result)
