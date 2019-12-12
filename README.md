# Ray Tracing + Anti-Aliasing
A ray tracer to handle opaque surfaces with lighting and shadows to realistically render a scene.

Change `SAMPLER_VALUE` on line #38 to alter anti-aliasing settings. (higher values take longer to render)\

Default is 3, which implies a 3x3 grid for direct supersampling.\

Setting it to a value of 1 turns AA off.

## Results

### Default Scene without Anti-Aliasing

![](raytracing-demo/images/default%20scene%20-%20no%20AA.png)

### Default Scene with Anti-Aliasing + `SAMPLER_VALUE 3`


![](raytracing-demo/images/default%20scene%20-%20sampler%20value%203%20-%20AA.png)


---
### Multiple Spheres Scene without Anti-Aliasing


![](raytracing-demo/images/spheres%20-%20no%20AA.png)


### Multiple Spheres Scene with Anti-Aliasing + `SAMPLER_VALUE 3`


![](raytracing-demo/images/spheres%20-%20sampler%20value%203%20-%20AA.png)


---
### Table Scene without Anti-Aliasing


![](raytracing-demo/images/table%20scene%20-%20no%20AA.png)


### Table Scene with Anti-Aliasing + `SAMPLER_VALUE 3`


![](raytracing-demo/images/table%20scene%20-%20sampler%20value%203%20-%20AA.png)


---
