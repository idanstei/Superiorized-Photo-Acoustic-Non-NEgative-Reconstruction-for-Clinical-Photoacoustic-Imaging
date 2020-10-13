# Superiorized-Photo-Acoustic-Non-NEgative-Reconstruction-for-Clinical-Photoacoustic-Imaging

This repository contains examples of the MATLAB code used to reconstruct photoacoustic images from the raw RF data using a fast superiorized conjugate gradient algorithm.
Results are published in our paper "Superiorized Photo-Acoustic Non-NEgative Reconstruction (SPANNER) for Clinical Photoacoustic Imaging".

The code is provided free for all to use. 
If you are publishing any work, where this code or variant of it has been used, please remember that it was obtained free of charge. 
We kindly ask that you reference our paper in the publication.

THE CODE IS PROVIDED "AS-IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OR CONDITIONS OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. 
IN NO EVENT SHALL THE AUTHORS AND/OR STANFORD UNIVERSITY BE LIABLE FOR ANY SPECIAL, INCIDENTAL, INDIRECT, OR CONSEQUENTIAL DAMAGES OF ANY KIND, OR DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA, OR PROFITS, WHETHER OR NOT THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES, AND/OR ON ANY THEORY OF LIABILITY ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS CODE.

MATLAB is a trademark of The MathWorks, Inc. Trademarks of other companies and/or organizations mentioned in this documentation and web-site appear for identification purposes only and are the property of their respective companies and/or organizations.

Contact:

Please contact Idan Steinberg PhD (idanstei@stanford.edu) for queries relating to this code.

Abstract:

Photoacoustic (PA) imaging can revolutionize medical ultrasound by augmenting it with molecular information. 
However, clinical translation of PA imaging remains a challenge due to the limited viewing angles and imaging depth. 
Described here is a new robust algorithm called Superiorized Photo-Acoustic Non-NEgative Reconstruction (SPANNER), designed to reconstruct PA images in real-time and to address these limitations. 
The method utilizes precise forward modeling of the PA propagation and reception of signals while accounting for the effects of acoustic absorption, element size, shape, and sensitivity, as well as the transducer's impulse response and directivity pattern. 
A fast superiorized conjugate gradient algorithm is used for inversion. 
SPANNER is compared to three reconstruction algorithms: delay-and-sum (DAS), universal back-projection (UBP), and model-based reconstruction (MBR). 
All four algorithms are applied to both simulations and experimental data acquired from tissue-mimicking phantoms, ex vivo tissue samples, and in vivo imaging of the prostates in patients. 
Simulations and phantom experiments highlight the ability of SPANNER to improve contrast to background ratio by up to 20 dB compared to all other algorithms, as well as a 3-fold increase in axial resolution compared to DAS and UBP. 
Applying SPANNER on contrast-enhanced PA images acquired from prostate cancer patients yielded a statistically significant difference before and after contrast agent administration, while the other three image reconstruction methods did not, thus highlighting SPANNER's performance in differentiating intrinsic from extrinsic PA signals and its ability to quantify PA signals from the contrast agent more accurately.

Reference:

Idan Steinberg, Jeesu Kim, Martin K. Schineider, Dongwoon Hyun, Aimen Zlitni, Sarah M. Hooper, Tal Klap, Geoffrey A. Sonn, Jeremy J. Dahl, Chulhong Kim, Senior Member and Sanjiv Sam Gambhir (2021) Superiorized Photo-Acoustic Non-NEgative Reconstruction (SPANNER) for Clinical Photoacoustic Imaging,....
