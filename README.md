![FIGURE1](https://github.com/user-attachments/assets/dea2a0c6-e449-43cb-bf4c-47bcd2a5a501)
# Stimuli and task
The task paradigm required maintenance of eye fixation within a 2° radius around the center of a monitor (Asus VG248QE: 24in, resolution 1920x1080, 144 Hz refresh rate). Upon ensuring the eye fixation, a stimulus is displayed in place of the fixation point for a 50-ms duration over a gray background, followed by a 600-ms blank interval. The stimuli were scaled to fit in a 7° by 7° square at the center of the monitor. Monkeys got rewarded for holding fixation every 1.5 seconds with droppings of sweet juice. The stimulus set contains 155 pictures from 4 different categories: face, body, artificial, and natural. 
Behavioral and recording method
We conducted experiments using acute extracellular recording in behaving rhesus macaques. Using two single tungsten electrodes (A-M Systems, FHC) simultaneously, we recorded the broadband neural activity of the prefrontal and inferior temporal cortices. Since the receptive field (RF) of IT neurons is relatively large and the RF of PFC neurons is not well-defined, we determined the spatial preference based on the PFC neuron’s RF. Before starting the main task, the RF of PFC neurons was determined by performing an RF-mapping task and monitoring the responses of the neurons to the spatial location of the stimulus. A reference grid, with holes 1 mm apart, was used inside all the recording chambers to guide electrode penetrations and localize them relative to structural MRI images. The neural response was amplified, quantized, and saved using a 40KHz data acquisition system (NikTek Systems).
# Neural data preprocessing
 Wide-band neural data was initially band-pass filtered between 300 to 3000 Hz to extract the high-frequency signals. By setting a multiplier of the median of the resulting waveform, the spikes were detected using the ROSS Toolbox 28. Using different electrodes for recording from different monkeys, we did not use the sorted spikes in the analysis and used the multiunit activity instead. Followed by spike detection, we chose a subset of recorded neurons using a selectivity criteria by which a unit was discarded if its post-stimulus average firing rate didn’t change significantly compared to its baseline activity for all of the four categories. The statistical significance was measured by a Wilcoxon signed rank test with a 5% confidence level.

#### Files Structure

```bash
SRC
    ├── Fano factor
                    ├── Category based
                    └── Stimuli based
    ├── GLM
    ├── Information theory
                        ├── Mutual Information
                        ├── Transfer Entropy
                        └── Factor Analysis
    ├── PSTH
    └── Machine learning methods
                            ├── LDA
                            ├── Generalization (Time-Time decoding)
                            └── SVM                          -       
```












