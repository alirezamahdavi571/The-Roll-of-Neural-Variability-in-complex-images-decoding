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
    ├── LFP
          └── PAC
    └── Machine learning methods
                            ├── LDA
                            ├── Generalization (Time-Time decoding)
                            └── SVM

                         -       
```

# **Transfer Entropy: Mathematical Explanation**

Transfer Entropy (TE) is a non-parametric statistical measure that quantifies the directional flow of information between two stochastic processes. It is grounded in information theory and can be particularly useful in understanding causal relationships in time series data.

## **Mathematical Definition**

Given two stochastic processes $X$ and $Y$, the Transfer Entropy from $X$ to $Y$, denoted as $T_{X \to Y}$, is defined as:

$$
T_{X \to Y} = \sum P(y_{t+1}, y_t^{(k)}, x_t^{(l)}) \log \frac{P(y_{t+1} \mid y_t^{(k)}, x_t^{(l)})}{P(y_{t+1} \mid y_t^{(k)})}
$$

Here:
- $y_{t+1}$: The state of process $Y$ at time $t+1$.
- $y_t^{(k)} = (y_t, y_{t-1}, \dots, y_{t-k+1})$: The past $k$ states of $Y$ (history of $Y$).
- $x_t^{(l)} = (x_t, x_{t-1}, \dots, x_{t-l+1})$: The past $l$ states of $X$ (history of $X$).
- $P(y_{t+1}, y_t^{(k)}, x_t^{(l)})$: The joint probability distribution of $y_{t+1}$, $y_t^{(k)}$, and $x_t^{(l)}$.
- $P(y_{t+1} \mid y_t^{(k)}, x_t^{(l)})$: The conditional probability of $y_{t+1}$ given $y_t^{(k)}$ and $x_t^{(l)}$.

## **Interpretation**
1. $T_{X \to Y}$ measures the reduction in uncertainty of $Y$'s future ($y_{t+1}$) by incorporating the history of $X$ ($x_t^{(l)}$), beyond what is already explained by the history of $Y$ itself ($y_t^{(k)}$).
2. A high $T_{X \to Y}$ value implies that $X$ has a strong influence on $Y$.

## **Key Features**
- **Directional**: Unlike correlation, TE quantifies the direction of information flow ($X \to Y$) and is not symmetric ($T_{X \to Y} \neq T_{Y \to X}$).
- **History Dependence**: By incorporating $k$- and $l$-length histories, TE captures lagged dependencies.
- **Non-linear Relationships**: TE can capture non-linear interactions between $X$ and $Y$.

## **Practical Applications**
- **Neuroscience**: Identifying causal interactions between brain regions.
- **Finance**: Understanding dependencies between market variables.
- **Engineering**: Analyzing complex systems such as power grids or ecological networks.

## **Simplified Example**
Assume two time series $X = \{x_1, x_2, ..., x_t\}$ and $Y = \{y_1, y_2, ..., y_t\}$. If the dynamics of $Y$ depend on both its own past and the past of $X$, TE helps determine how much of $Y$'s future behavior is influenced by $X$'s past.



# Normalized Direct Phase-Amplitude Coupling (ndPAC) Method

The **Normalized Direct Phase-Amplitude Coupling (ndPAC)** method, introduced by Ozkurt et al. (2012), quantifies the interaction between the phase of low-frequency oscillations and the amplitude of high-frequency oscillations. It is widely applied in neuroscience to study neural rhythm interactions.

---

## How ndPAC Works

### 1. Signal Filtering
The input signal is bandpass-filtered to isolate:
- **Low-frequency component** (\( f_\text{phase} \))
- **High-frequency component** (\( f_\text{amplitude} \))

### 2. Hilbert Transform
The analytic signals for each filtered component are computed using the Hilbert transform:

\[
z_\text{phase}(t) = a_\text{phase}(t) e^{i\phi(t)}
\]

\[
z_\text{amp}(t) = a_\text{amp}(t) e^{i\theta(t)}
\]

Where:
- \( \phi(t) \): Phase of the low-frequency component.
- \( a_\text{amp}(t) \): Amplitude envelope of the high-frequency component.

### 3. Phase-Amplitude Coupling (PAC) Calculation
The PAC value is computed as:

\[
\text{PAC} = \left| \frac{1}{N} \sum_{t=1}^N a_\text{amp}(t) e^{i\phi(t)} \right|
\]

Here:
- \( a_\text{amp}(t) \): High-frequency amplitude envelope modulated by the low-frequency phase \( \phi(t) \).
- \( N \): Number of time samples.

### 4. Normalization
The PAC is normalized to reduce bias in amplitude modulation estimates. This is achieved either by comparing PAC against surrogate data or through analytical normalization methods.

---

## Key Features of ndPAC

- **Direct Calculation:** Bypasses the need for binning or averaging, improving computational efficiency.
- **Normalization:** Provides unbiased and reliable coupling estimates.

---

This method enables precise and normalized quantification of phase-amplitude coupling, making it an essential tool for analyzing neural dynamics and studying brain rhythms.









