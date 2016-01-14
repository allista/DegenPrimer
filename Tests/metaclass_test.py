'''
Created on Mar 5, 2014

@author: allis
'''

class Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

class MyClass(object):
    __metaclass__ = Singleton
    
    
o = MyClass()