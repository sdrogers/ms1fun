"""
Tests for some of the trickier parts of the grouping algorithm.
"""

def singleton_test(groups):
    """
    Checks whether every singleton group has an M+H associated with them. 
    They should.
    :param groups:  list of peak groups that need to be checked.
    :return:        True:   all singletons either have M+H, or there are no singleton groups. 
                    False:  they don't.
    """
    singletons = 0
    for group in groups:
        if len(group.members) == 1:
            transformation = group.members[0][1]  # if one member - always [0]
            singletons += 1
            if transformation.name != 'M+H':
                return False
    
    if singletons < 1:
        raise UserWarning("There were no singletons.")

    return True
