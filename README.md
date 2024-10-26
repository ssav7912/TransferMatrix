# Realtime Transfer Matrix Layered Materials 
An extension of MiniEngine by Team Minigraph at Microsoft
Based on research in the following paper: 
Joël Randrianandrasana, Patrick Callet, and Laurent Lucas. 2021. Transfer matrix based layered materials rendering. ACM Trans. Graph. 40, 4, Article 177 (August 2021), 16 pages. https://doi.org/10.1145/3450626.3459859

  

## Getting started:
* Open ModelViewer/ModelViewer_VS16.sln
* Select configuration: Debug (full validation), Profile (instrumented), Release
* Select platform
* Build and run

## Running the TransferMatrix shading model:
* Make sure the 'Required Resources' .zip from the github releases has been extracted to the ModelViewer directory.
* Add the `-model TestSphere.glb -TransferMatrix 1` flags to the commandline when starting. The geometry may be swapped for any other GLTF in the ModelViewer directory if desired.

## Controls:
* forward/backward/strafe: left thumbstick or WASD (FPS controls)
* up/down: triggers or E/Q
* yaw/pitch: right thumbstick or mouse
* toggle slow movement: click left thumbstick or lshift
* open debug menu: back button or backspace
* navigate debug menu: dpad or arrow keys
* toggle debug menu item: A button or return
* adjust debug menu value: dpad left/right or left/right arrow keys

## glTF 2.0 Support:

Get sample assets from https://github.com/KhronosGroup/glTF-Sample-Models or make your own.

* Place asset folder underneath the ModelViewer folder
* Add to command line "-model [relative path to gltf or glb file]"
* Example:  ModelViewer.exe -model SciFiHelmet/glTF/SciFiHelmet.gltf

Notes:  Some IBL cube maps are provided for physically based rendering in the ModelViewer/Textures folder.  ModelViewer automatically detects cube maps in this folder so you can add your own.  They must follow the same naming convention for diffuse and specular maps.  You can change the active environment map in the tweak menu under the "ModelViewer" category.

DirectXMesh and DirectXTex are used for compiling content the first time it is loaded.  For example, image files will have mip maps generated, encoded to a block compressed format, and then saved as a DDS for subsequent loading.  This means that first-time loading takes a little longer.


## Development Notes:
A number of preprocessor macros are available for enabling/disabling various features, in order to reproduce experiments. These are:
```
USE_KARIS_FGD=0 //if 1, uses the Karis Split-Sum LUT. Mutually Exclusive with USE_BELCOUR_FGD
USE_BELCOUR_FGD=1 //if 1, uses the Belcour Split-Sum LUT. Mutually Exclusive with USE_KARIS_FGD. 
//if both above are 0, the original resampled FGD LUT is used.


USEFP16=0 //Enables half precision types.
USE_D_KARIS=0 //If 1, Use D term from 'Real Shading in Unreal Engine 4' (Karis, 2013)

USE_EARL_G2=0 //If 1, Use Earl G2 Approximation. Mutually exclusive with USE_SMITH_G2
USE_SMITH_G2=0 //If 1, Use smith G2 function.
SCHLICK_G=0 //If 1, Use Schlick G2 function.
//if all three above are 0, uses the height-correlated geometry term from 'Understanding The Masking Shadowing Function in Microfacet Based BRDFs' (Heitz, 2014)
DISABLE_TIR=0 //If 1, disables TIR corrections.
ANALYTIC_TIR=0 //If 1, uses the analytic TIR approximation developed in this research.
```