def Settings( **kwargs ):
  return {
    'flags': [ '-x', 'c++', '-Wall', '-Wextra', '-Werror', '-include', '/usr/local/MATLAB/R2021a/extern/include/mex.h',
        '-include', '/usr/local/MATLAB/R2021a/extern/include/mexAdapter.hpp'],
  }
