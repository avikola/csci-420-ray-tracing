Avishkar Kolahalu

------------------


Input files available to test:
test2.txt
spheres.txt
table.txt



Feature:                                 Status: finish? (yes/no)

-------------------------------------    -------------------------

1) Ray tracing triangles                  yes


2) Ray tracing sphere                     yes


3) Triangle Phong Shading                 yes


4) Sphere Phong Shading                   yes


5) Shadows rays                           yes


6) Still images                           yes


Image Details:
(images located in the 'screenshots' folder)

000.jpg - test2.txt
001.jpg - test2.txt - with anti-aliasing

002.jpg - spheres.txt
003.jpg - spheres.txt - with anti-aliasing

004.jpg - table.txt
005.jpg - table.txt - with anti-aliasing


• Implemented Anti-Aliasing:
	• Used direct supersampling method to deal with aliasing.
	• Current grid size is 3x3.
	• Change SAMPLER_VALUE to adjust the grid size.
		Note: higher values take a longer time.
