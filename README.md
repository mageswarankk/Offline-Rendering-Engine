## 3D Path-Tracing Rendering Engine using Dielectric Microfacet BSDF (with GGX Normal Distribution)

This project was my Final Project for CS190I: Introduction to Offline Rendering taught by Professor Lingqi Yan at UCSB

For this project, I created a rendered image by employing path tracing enhanced with multiple importance sampling. The shading effect was achieved through microfacet models, utilizing both the standard microfacet Bidirectional Reflectance Distribution Function (BRDF) for reflective materials, such as iron and gold, and a Bidirectional Scattering Distribution Function (BSDF) to capture both reflectance and transmittance properties. This approach allowed for the realistic simulation of non-conductive materials like glass, air, water, and ice through the implementation of a dielectric BSDF in my rendering engine. This dielectric BSDF combines a BRDF for reflectance and a BTDF (Bidirectional Transmittance Distribution Function) for transmittance, with the interplay between these components governed by the Fresnel Term as determined by Snell’s Law. For microfacet direction sampling, I opted for the GGX distribution, which offers sharper peaks and extended tails, thereby enabling a more nuanced representation of material roughness.

In developing the code for this project, I built upon the framework provided in Assignment 4, which was based on the Muni Rendering Toolchain in C++. Starting with my previous work that included an implementation of the GGX distribution as a bonus, I streamlined the code by removing redundant functions, specifically those associated with light sampling and the Beckmann Normal Distribution that were not utilized in this project. Additionally, I introduced a new struct, $\texttt{Dielectric}$, designed to encapsulate the material-specific functions necessary for a comprehensive set of $\texttt{eval()}$, $\texttt{pdf()}$, and $\texttt{sample()}$ functions, further enhancing the renderer's ability to simulate complex material interactions accurately.

Finally here is my final rendered image, produced in 1080 x 1080 resolution with 128 samples per pixel. It contains two bunny images made of a Dielectric material simulating glass (with index of refraction 1.55) with the smaller bunny being of considerably higher roughness and also tinted in a shade of blue-green using the roughness and attenuation parameter. (Note the colors may be misleading since the light source is not white color - a stylistic choice I made). The walls are produced with lambertian material, and the floor is an iron microfacet material. The walls are derived from the Cornell Box.

Raw Image:
![raw_image](https://github.com/mageswarankk/offline-rendering-engine/assets/73411210/5e9424a9-1677-4a9b-9328-36381dccc865)

Denoised Image:
![denoised_image](https://github.com/mageswarankk/offline-rendering-engine/assets/73411210/e085abaf-6f8c-459c-86c0-363db04ed4b7)
