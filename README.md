# PichiaCLM - PyTorch Implementation

This project is the PyTorch port and optimized version of **PichiaCLM** (a deep learning-based codon optimization model for *Pichia pastoris*).

> **🌟 Original Source**
> The core architecture and original data processing logic of this project are derived from [NHarini-1995/PichiaCLM](https://github.com/NHarini-1995/PichiaCLM). The original code was written in Keras/TensorFlow; this project has completely rewritten it in PyTorch, featuring deep mathematical and biological rule optimizations in the inference logic.

---

## 🛠️ Modifications & Improvements

During the migration from Keras to PyTorch, we not only replicated the original Multi-Task Seq2Seq architecture but also implemented specific fixes and upgrades for issues exposed during actual inference:

### 1. Complete Framework Migration
* **Model Construction**: Replaced the original Keras layers with PyTorch's `nn.Embedding` and `nn.GRU`, and manually implemented a **Dot-Product Attention** mechanism strictly consistent with the original version.
* **Loss Function**: Utilized `nn.CrossEntropyLoss` with the `ignore_index=0` parameter to natively and efficiently achieve the sequence padding masking effect that was handled by `mask_zero=True` in Keras.
* **Training Loop**: Built a flexible, custom PyTorch training loop and incorporated gradient clipping (`clip_grad_norm_`), NaN/Inf loss interception, and a dynamic learning rate scheduler based on validation performance (`ReduceLROnPlateau`).

### 2. Fixing "Prefix Hallucination" and Alignment Shift
In early Seq2Seq tests, we discovered that the model would "hallucinate" when faced with sequences that did not start with an initiation codon (e.g., `M`), swallowing the first few amino acids and outputting entirely irrelevant codon combinations.
* **Refactored to Padding-Free Decoding**: Eliminated fixed-length Padding (zero-padding) during the inference phase, ensuring that the hidden states of the bidirectional GRU are not polluted by the noise of long sequences of zeros.
* **Pure Step-by-Step Autoregressive Decoding**: Changed the original approach of feeding the entire sequence at once to a strict, standard machine translation inference paradigm—feeding only a single Token at a time and updating the Hidden State. This allows the Attention mechanism to focus intensely on the true sequence length.

### 3. Hard Biological Rule Masking — Core Breakthrough 🚀
Standard language models carry the risk of probability out-of-bounds errors when translating amino acids. To achieve **100% biological translation fidelity**, we introduced a strict rule mask:
* **Dictionary Mapping Constraint**: When the decoder outputs the probabilities (Logits) for the next codon, the code intercepts the amino acid ID corresponding to the current position in real-time and queries the set of valid codons.
* **`-inf` Masking Mechanism**: For all invalid codons in the vocabulary that do not belong to that amino acid, the model forcefully adds negative infinity (`-inf`) to their Logits.
* **Effect**: This mechanism forces the model to **only select the codon deemed most suitable for the host's (*Pichia*) preferences by the Attention mechanism, chosen exclusively from the valid codon options**. This completely eradicates alignment issues like skipping, missing, or mistranslating, achieving a perfect integration of AI prediction and biological rules.

---

## 🚀 Quick Start

### Dependencies
* Python 3.8+
* PyTorch 2.0+
* Pandas, NumPy

### Core Files Description
* `/Model_PichiaCLM/Training/AllData/`: Training data preparation.
* `Train_Test_Splitting-Torch.ipynb`: Data cleaning, unified vocabulary construction (SOS, EOS, UNK, PAD), and Tensor serialization code.
* `Training-Torch.ipynb`: Contains the PyTorch model architecture, Multi-task Training loop, and the final optimized inference script with Hard Rule Masking.
* `./Model_PichiaCLM/Training/PichiaData/2Target_AllData/Arch1-0404.weights.pt`: Trained weights. Val Accuracy: CDS = 88.39% | AA = 99.99%.

### Inference Demo Snippet
The model can accurately generate codons based on the amino acid sequence, incorporating the expression preferences of *Pichia pastoris*:

```python
# Initialize and load the model
model = MultiTaskSeq2Seq().to(device)
model.load_state_dict(torch.load("Arch1-0404.weights.pt"))

# Input target amino acid sequence
input_aa_string = "DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADES"

# Execute decoding with biological constraints
predicted_indices = translate_aa_to_cds(model, aa_tensor, aa_id_to_codon_ids, device)
```
### Additional Notes
The CDS sequence inferred using the trained weights may not exactly match the data presented in the original paper.

## Acknowledgments
Thanks to the original author NHarini-1995 for the excellent model design and the Pichia pastoris dataset. This PyTorch version aims to provide the community with a modernized deep learning implementation that is easier to debug, extend, and develop further.

Coded with Gemini 3.1 Pro - 20260406