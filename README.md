# Ray Tracing + Anti-Aliasing
A ray tracer to handle opaque surfaces with lighting and shadows to realistically render a scene.

Change `SAMPLER_VALUE` on line #38 to alter anti-aliasing settings.\
Setting it to a value of 1 turns it off.

## Results

### Default Scene without Anti-Aliasing

![](raytracing-demo/images/default%20scene%20-%20no%20AA.png)

### Default Scene with Anti-Aliasing + `SAMPLER_VALUE 3`


![](raytracing-demo/images/default%20scene%20-%20qsampler%20value%203%20-%20AA.png)


### Multiple Spheres Scene without Anti-Aliasing


![](raytracing-demo/images/spheres%20-%20no%20AA.png)


### Multiple Spheres Scene with Anti-Aliasing + `SAMPLER_VALUE 3`


![](raytracing-demo/images/spheres%20-%20qsampler%20value%203%20-%20AA.png)


### Table Scene without Anti-Aliasing


![](raytracing-demo/images/table%20scene%20-%20no%20AA.png)


### Table Scene with Anti-Aliasing + `SAMPLER_VALUE 3`


![](raytracing-demo/images/table%20scene%20-%20sampler%20value%203%20-%20AA.png)

