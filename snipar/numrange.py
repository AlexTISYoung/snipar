import argparse, re
def parseNumRange(string):
    """reads either a int or a range"""
    match_range = re.match(r' *(\d+) *- *(\d+) *', string)
    match_list = re.match(r' *(\d+) *', string)
    if match_range:
        start = int(match_range.group(1))
        end = int(match_range.group(2))
        result = [str(n) for n in range(start, end+1)]
    elif match_list:
        result = match_list.group(0)
    else:
        raise Exception(f"{string} is neither a range of the form x-y nor a list of integers of the form x y z")
    return result

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