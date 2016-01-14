'''
Created on Mar 8, 2014

@author: allis
'''


def countable(func):
    def count_func(*args, **kwargs):
        func.func_globals['counter'] = kwargs.pop('counter', 1) 
        return func(*args, **kwargs)
    return count_func


if __name__ == '__main__':
    @countable
    def f(a, b):
        global counter
        counter += 1
        return a+b+counter
    
    
    print f(*(1,2), counter=2)