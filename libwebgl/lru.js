/* Least Recently Used (LRU) cache management at client side */
class LRU 
{
  /*
  .. create an LRU object with 'max' as the limit of cache.
  .. Create a hash map for insert/probing. when map size reaches
  .. 'max' it evicts
  */
  constructor(max)
  {
    this.max = max;
    this.map = new Map();
  }

     
  get(k)
  {
    if (!this.map.has(k))
      return undefined;
    const v = this.map.get(k);
    /* "recently" used */
    this.map.delete(k); this.map.set(k, v);
    return v;
  }

  set(k, v, onEvict)
  {
    this.map.set(k, v);
    if (this.map.size > this.max) 
    {
      const oldest = this.map.keys().next().value;
      onEvict?.(oldest, this.map.get(oldest));
      this.map.delete(oldest);
    }
  }
}
