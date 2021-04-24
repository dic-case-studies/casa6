
def create_error_string( errors ):
    variable = None
    value = None
    vec = False
    if type(errors) == dict:
        for k,v in errors.items():
            if type(v) == list:
                if v[0] == 'path not found':
                    value = k
                if v[0] == 'path vector element not found':
                    value = k
                    vec = True
                if v[0] == 'must be of cReqPath type':
                    variable = k
    if variable is not None and value is not None:
        return ( "all elements of the %s vector parameter must be a path that exists (at least one element of \"%s\" does not exist)"
                 if vec else "the %s parameter must be a path that exists ('%s' does not exist)" ) % (variable, value)
    else:
        return str(errors)
