# Superiorized-Photo-Acoustic-Non-NEgative-Reconstruction-for-Clinical-Photoacoustic-Imaging
Photoacoustic (PA) imaging can revolutionize medical ultrasound by augmenting it with molecular information. However, clinical translation of PA imaging remains a challenge due to the limited viewing angles and imaging depth. Described here is a new robust algorithm called Superiorized Photo-Acoustic Non-NEgative Reconstruction (SPANNER), designed to reconstruct PA images in real-time and to address these limitations. The method utilizes precise forward modeling of the PA propagation and reception of signals while accounting for the effects of acoustic absorption, element size, shape, and sensitivity, as well as the transducer's impulse response and directivity pattern. A fast superiorized conjugate gradient algorithm is used for inversion. SPANNER is compared to three reconstruction algorithms: delay-and-sum (DAS), universal back-projection (UBP), and model-based reconstruction (MBR). All four algorithms are applied to both simulations and experimental data acquired from tissue-mimicking phantoms, ex vivo tissue samples, and in vivo imaging of the prostates in patients. Simulations and phantom experiments highlight the ability of SPANNER to improve contrast to background ratio by up to 20 dB compared to all other algorithms, as well as a 3-fold increase in axial resolution compared to DAS and UBP. Applying SPANNER on contrast-enhanced PA images acquired from prostate cancer patients yielded a statistically significant difference before and after contrast agent administration, while the other three image reconstruction methods did not, thus highlighting SPANNER's performance in differentiating intrinsic from extrinsic PA signals and its ability to quantify PA signals from the contrast agent more accurately.
