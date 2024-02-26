import functools

# https://realpython.com/primer-on-python-decorators/#caching-expensive-functions

def cache(func):
    """Keep a cache of previous function calls"""
    @functools.wraps(func)
    def wrapper_cache(*args, **kwargs):
        cache_key = args + tuple(kwargs.items())
        if cache_key not in wrapper_cache.cache:
            print("cache miss")
            wrapper_cache.cache[cache_key] = func(*args, **kwargs)
        else:
            print("cache hit")
        return wrapper_cache.cache[cache_key]
    wrapper_cache.cache = dict()
    return wrapper_cache