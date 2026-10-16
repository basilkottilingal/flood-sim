# client side

This is the front end that uses webgl to render geotiff files by communicating with the wsserver.

## Still In Progress...

1. Need to decide on how to split work b/w client & server
   [Since you need server to do a lot of computation, you can employ the server to handle
   a lot of webgl/webgpu computation]
2. Need to convert 2D [[lat,lon]=> alt] raster to [[x,y,z]] or [[x,y] => z] raster
3. Camera, zoom, ..
4. update Tiles in the scene for a given zoom. May use local LRU cache and if not available
   may ask the server
  4.1 tile management [LRU + request]
  4.2 tile decoding [MVT ??]
5. style spec. Hillshading, layers (vegetation), etc
6. Fill, Line, 
7. Symbols (name, tags, ..)
8. Options and fetch

## Read on

Webgpu rendering, falling back to webgl 2.0. (Note that webgl 1.0 cannot support three.js. However
webgl 2.0 is universal these days, it is very rare to see a system/browser that only has wegbl 1.0)

three.js, webgl, webgpu, 

## Designing

As I am new to webgl, ...


