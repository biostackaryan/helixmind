# Helix Mind – Bioinformatics Toolkit

Helix Mind is an easy-to-use bioinformatics web app that brings together multiple tools for sequence analysis, literature search, and 3D visualization — all in one place.

## 🔬 Features
- **FASTA File Analysis** – Upload your FASTA file to view sequence statistics like GC content and length.
- **BLAST Search** – Run BLAST directly in the app to find sequence similarities.
- **PubMed Search** – Look up relevant scientific articles quickly.
- **KEGG Pathways** – Explore metabolic and molecular pathways.
- **3D Structure Viewer** – Visualize protein and nucleic acid structures by PDB ID.
- **AI Assistant** – Get intelligent responses powered by ChatGPT (if API key is provided).

## 🖥 How to Use
1. Upload your FASTA file or enter sequence data.
2. Select the tool you want to use (BLAST, PubMed, KEGG, etc.).
3. View results instantly in the browser.
4. For AI Assistant, enter your OpenAI API key in the input field.

## 📦 Installation (Local)
If you want to run it locally instead of Hugging Face:
```bash
git clone https://github.com/biostackaryan/helixmind.git
cd helixmind
pip install -r requirements.txt
python app.py
