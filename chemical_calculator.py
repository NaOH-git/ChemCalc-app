import streamlit as st
import requests

st.set_page_config(page_title="Chemical Calculator", page_icon="⚗️")
st.title("⚗️ Chemical Calculator with Unit Conversion")

# --- Session State for Persistent Outputs ---
if 'mass_result' not in st.session_state:
    st.session_state.mass_result = None
if 'gl_result' not in st.session_state:
    st.session_state.gl_result = None
if 'molecular_weight' not in st.session_state:
    st.session_state.molecular_weight = 0.0

# --- Helper functions ---
def get_multiplier(unit):
    return {
        "L": 1,
        "mL": 1e-3,
        "µL": 1e-6,
        "M": 1,
        "mM": 1e-3,
        "µM": 1e-6,
        "nM": 1e-9,
        "g": 1,
        "mg": 1e-3,
        "µg": 1e-6,
        "ng": 1e-9
    }[unit]

def format_mass(value_in_g):
    units = ["g", "mg", "µg", "ng"]
    factors = [1, 1e-3, 1e-6, 1e-9]
    for unit, factor in zip(units, factors):
        if value_in_g >= factor:
            return f"{value_in_g / factor:.4f} {unit}"
    return f"{value_in_g / 1e-9:.4f} ng"

def fetch_molecular_weight(name):
    try:
        # First: get CID
        cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
        cid_response = requests.get(cid_url)
        cid_response.raise_for_status()
        cid_data = cid_response.json()
        cid = cid_data['IdentifierList']['CID'][0]

        # Second: get MolecularWeight by CID
        mw_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight/JSON"
        mw_response = requests.get(mw_url)
        mw_response.raise_for_status()
        mw_data = mw_response.json()
        mw = mw_data['PropertyTable']['Properties'][0]['MolecularWeight']
        return float(mw), cid
    except:
        return None, None

def fetch_uniprot_mw(name):
    try:
        search_url = f"https://rest.uniprot.org/uniprotkb/search?query={name}&fields=accession,mass,protein_name&format=json&size=1"
        response = requests.get(search_url)
        response.raise_for_status()
        data = response.json()
        entry = data['results'][0]
        mass = entry['proteinDescription']['recommendedName']['fullName']['value']
        mw = entry['sequence']['mass']  # in Da
        accession = entry['primaryAccession']
        return mw, accession
    except:
        return None, None

# --- Molecular Weight Lookup ---
st.subheader("Molecular Weight Lookup")
compound_name = st.text_input("Compound Name (e.g., glucose, NaCl.., (english writing))")
search_as_protein = st.checkbox("Search as protein (UniProt)")

if st.button("Lookup Molecular Weight"):
    if search_as_protein:
        mw_result, acc = fetch_uniprot_mw(compound_name)
        if mw_result:
            mw_kDa = mw_result / 1000
            st.success(f"Estimated Mass of {compound_name}: {mw_kDa:.2f} kDa")
            st.session_state.molecular_weight = mw_result
            st.markdown(f"[🔗 View on UniProt](https://www.uniprot.org/uniprotkb/{acc})")
        else:
            st.error("Protein not found on UniProt.")
    else:
        mw_result, cid_result = fetch_molecular_weight(compound_name)
        if mw_result:
            st.success(f"Molecular Weight of {compound_name}: {mw_result:.2f} g/mol")
            st.session_state.molecular_weight = mw_result
            st.markdown(f"[🔗 View on PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/{cid_result})")
        else:
            st.error("Compound not found or error retrieving data.")

# --- Mass Calculator ---
st.subheader("Mass Calculator")
col1, col2 = st.columns([3, 1])
with col1:
    mw = st.number_input("Molecular Weight", min_value=0.0, format="%f", value=st.session_state.molecular_weight)
    st.session_state.molecular_weight = mw
with col2:
    st.markdown("**g/mol**")

col3, col4 = st.columns([3, 1])
with col3:
    volume = st.number_input("Volume", min_value=0.0, format="%f")
with col4:
    volume_unit = st.selectbox("Volume Unit", ["L", "mL", "µL"], key="vol_unit")

col5, col6 = st.columns([3, 1])
with col5:
    concentration = st.number_input("Concentration", min_value=0.0, format="%f")
with col6:
    concentration_unit = st.selectbox("Concentration Unit", ["M", "mM", "µM", "nM"], key="conc_unit")

if st.button("Calculate Required Mass"):
    try:
        vol_L = volume * get_multiplier(volume_unit)
        conc_M = concentration * get_multiplier(concentration_unit)
        mass_g = mw * conc_M * vol_L
        breakdown = f"{concentration} {concentration_unit} * {volume} {volume_unit} * {mw} g/mol"
        st.session_state.mass_result = f"Mass required: {format_mass(mass_g)}"
        st.session_state.mass_breakdown = breakdown
    except Exception as e:
        st.session_state.mass_result = f"Error: {e}"

# --- Molarity to g/L Conversion ---
st.subheader("Convert Molarity to g/L")
if st.button("Convert to g/L"):
    try:
        conc_M = concentration * get_multiplier(concentration_unit)
        gl = mw * conc_M
        st.session_state.gl_result = f"Equivalent concentration: {gl:.4f} g/L"
    except Exception as e:
        st.session_state.gl_result = f"Error: {e}"

# --- Display Combined Results ---
if st.session_state.mass_result:
    st.success(st.session_state.mass_result)
    if 'mass_breakdown' in st.session_state:
        st.caption(f"Calculation: {st.session_state.mass_breakdown}")

if st.session_state.gl_result:
    st.info(st.session_state.gl_result)

# --- Reverse Calculation ---
st.subheader("Calculate Molarity from Mass and Volume")
col7, col8 = st.columns([3, 1])
with col7:
    mass_input = st.number_input("Mass", min_value=0.0, format="%f")
with col8:
    mass_unit = st.selectbox("", ["g", "mg", "µg", "ng"])

col9, col10 = st.columns([3, 1])
with col9:
    volume2 = st.number_input("Volume (for reverse)", min_value=0.0, format="%f")
with col10:
    volume_unit2 = st.selectbox(" ", ["L", "mL", "µL"], key="vol_unit2")

if st.button("Calculate Molarity"):
    try:
        vol_L = volume2 * get_multiplier(volume_unit2)
        mass_g = mass_input * get_multiplier(mass_unit)
        molarity = mass_g / mw / vol_L
        st.success(f"Molarity: {molarity:.6f} mol/L")
    except Exception as e:
        st.error(f"Error: {e}")
