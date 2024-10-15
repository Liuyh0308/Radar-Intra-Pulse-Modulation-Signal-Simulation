# Radar-Intra-Pulse-Modulation-Signal-Simulation
Code for Author's Article：The Research of Intra-Pulse Modulated Signal Recognition of Radar Emitter under Few-Shot Learning Condition Based on Multimodal Fusion.
url：https://doi.org/10.3390/electronics13204045

Mentions：
1. If you want to get the radar in-pulse modulation signal simulation dataset, run the cdwtfi_data_generation code and use the WVD time-frequency transform to generate a 2D time-frequency plot, and at the same time save the signal sequence of the radar in-pulse modulation signal as a mat format file as well。
2. If you just want to generate a sequence of radar signals and save them, run mat_data_gen. mat_data_gen_5pathdelay is a signal simulation environment that incorporates time delays and other factors.
3. Make sure you download the “tftb” signal time-frequency analysis toolkit and introduce the path to your MATLAB before running the code.
4. The download link for the tftb toolkit is: http://www.nongnu.org/tftb/
5. Some places may need to be commented or un-commented to achieve your needs, such as getting the time-frequency plot dataset and saving it automatically, for details, please refer to the comments of the specific code and instructions for use.
