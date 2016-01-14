# Copyright (C) 2012 Allis Tauri <allista@gmail.com>
# 
# degen_primer is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# degen_primer is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
Created on Feb 26, 2014

@author: Allis Tauri <allista@gmail.com>
'''

import os, atexit
from shelve import DbfilenameShelf
from tempfile import mkstemp, gettempdir
from UMP import UManager 

class _SetManager(UManager): pass
_SetManager.register('set', set)
_set_manager = _SetManager()
_set_manager.start()

_tmp_files = _set_manager.set()

def register_tmp_file(path):
    global _tmp_files
    _tmp_files.add(path)
#end def

def unregister_tmp_file(path):
    global _tmp_files
    try: _tmp_files.remove(path)
    except: pass
#end def

def cleanup_file(fpath):
    try: os.unlink(fpath)
    except:
        #this means a dbm module was used by anydbm (in DbfilenameShelf
        #that makes several files with the base self.filename and different
        #extensions
        tmp_dir = os.path.dirname(fpath)
        fname   = os.path.basename(fpath)
        for entry in os.listdir(tmp_dir):
            if entry.startswith(fname):
                try: os.unlink(entry)
                except: pass
#end def

@atexit.register
def clean_tmp_files():
    global _tmp_files
    for fpath in _tmp_files._getvalue(): cleanup_file(fpath)
#end def


class tmpDict(DbfilenameShelf):
    def __init__(self, fname=None, tmpdir=None, persistent=False):
        if fname is None:
            #get tmp directory
            if tmpdir is None:
                self._tmpdir = gettempdir()
            else: self._tmpdir = tmpdir
            #create tmp file for db
            fd, self.filename = mkstemp('', 'DP-PCR-', dir=self._tmpdir)
            os.close(fd); os.unlink(self.filename)
            #create a shelf in the db
            DbfilenameShelf.__init__(self, self.filename, flag='n', protocol=-1)
            if not persistent: register_tmp_file(self.filename)
        else:
            self._tmpdir  = os.path.dirname(fname)
            self.filename = fname
            DbfilenameShelf.__init__(self, self.filename, flag='w', protocol=-1)
    #end def
    
    
    def cleanup(self):
        try: self.close()
        except: pass
        cleanup_file(self.filename)
    #end def
#end def


class roDict(DbfilenameShelf):
    def __init__(self, filename):
        self.filename = filename
        DbfilenameShelf.__init__(self, filename, flag='r', protocol=-1)
        
    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()
#end class


from collections import Sequence
class tupleView(Sequence):
    def __init__(self, keys, db):
        self._db = db
        self._keys = keys
        
    def __len__(self):
        return len(self._keys)
        
    def __getitem__(self, index):
        key = self._keys[index]
        if isinstance(key, str):
            return self._db[key]
        elif isinstance(key, Sequence):
            return tuple(self._db[k] for k in key)
        raise ValueError('Unsupported key type: %s' % type(key)) 
    #end def
    
    def __reduce__(self):
        return _unpickle_tupleView, (self._keys, self._db.filename)
#end class 

def _unpickle_tupleView(keys, db_path):
    return tupleView(keys, roDict(db_path))


def shelf_result(func, *args, **kwargs):
    '''Decorator that saves result of func(*args) into a temporary shelf 
    under "result" key and returns the path to it's database'''
    def s_func(*args, **kwargs):
        result = func(*args, **kwargs)
        if result is None: return None
        d = tmpDict(persistent=True)
        d['result'] = result 
        d.close()
        return d.filename
    s_func.__name__ = func.__name__+'_to_shelf'
    return s_func
#end def



#tests
if __name__ == '__main__':
    from time import sleep
    import cPickle as pickle
    
    d = tmpDict()
    d['a'] = (1,2,3,4,5)
    d['b'] = dict(c=(6,7), d=(8,9))
    d['c'] = 'asdfas'
    d['d'] = 1245134613
    fname = d.filename
    print d
    print fname
    print 'File exists:', os.path.isfile(fname)
    tv = tupleView(('a','b','c','d'), d)
    tv1 = pickle.loads(pickle.dumps(tv, protocol=-1))
    print tv[0], tv[1]
    print tv[2:4]
    del d
    print 'Was file deleted?', (not os.path.isfile(fname))