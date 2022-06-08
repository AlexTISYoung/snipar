def get_parser_doc(parser):
    doc = ""
    for action in parser._actions:
        options = str(action.option_strings)[1:-1]
        default = action.default
        type = ""
        if action.type:
            if "'" in str(action.type):
                type = str(action.type).split("'")[1]
        
        help = action.help
        default_substring = ""
        if default:
            default_substring = f", default={default}"
        
        arg_doc = f"""    {options} : {type}{default_substring}
            {help}

"""
        doc += arg_doc
    
    return doc