def Settings( **kwargs ):
  return {
    'flags': [ '-x', 'cuda', '--cuda-gpu-arch=sm_50',
        '-L/usr/local/cuda/lib64', '-lcudart_static', '-ldl', '-lrt', '-pthread'],
  }
