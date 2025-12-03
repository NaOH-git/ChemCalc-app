# chemical_calculator_v2_2_1.py
import streamlit as st
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import isfinite

st.set_page_config(page_title="Chemical Calculator", page_icon="⚗️", layout="wide")

# -------------------------
# Helper functions
# -------------------------
def get_multiplier(unit):
    return {
        "L": 1.0, "mL": 1e-3, "µL": 1e-6,
        "M": 1.0, "mM": 1e-3, "µM": 1e-6, "nM": 1e-9,
        "mg/mL": 1.0, "µg/mL": 1e-3, "ng/mL": 1e-6
    }[unit]

def is_mass_unit(unit):
    return unit in ("mg/mL", "µg/mL", "ng/mL")

def is_molar_unit(unit):
    return unit in ("M", "mM", "µM", "nM")

def conc_to_g_per_L(value, unit, mw=None):
    """
    Convert concentration value (in given unit) to g/L.
    - For mass-based units (mg/mL, µg/mL, ng/mL) convert directly.
      Note: 1 mg/mL == 1 g/L.
    - For molar units, MW (g/mol) is required: (mol/L) * (g/mol) = g/L
    """
    if is_mass_unit(unit):
        # convert mass conc to g/L using factors stored in get_multiplier:
        # value * factor => value in (g/L) because mg/mL -> 1 g/L etc.
        factor = get_multiplier(unit)
        return float(value) * factor
    elif is_molar_unit(unit):
        if mw is None or mw <= 0:
            raise ValueError("Molecular weight required for molar → mass conversions.")
        # get mol/L
        mol_per_L = float(value) * get_multiplier(unit)
        return mol_per_L * float(mw)
    else:
        raise ValueError(f"Unknown concentration unit: {unit}")

def format_mass(value_in_g):
    units = ["g", "mg", "µg", "ng"]
    factors = [1, 1e-3, 1e-6, 1e-9]
    for u, f in zip(units, factors):
        if value_in_g >= f:
            return f"{value_in_g / f:.4f} {u}"
    return f"{value_in_g / 1e-9:.4f} ng"

def fetch_molecular_weight_and_name(name):
    """
    Returns tuple: (mw_float_or_None, cid_or_None, preferred_name_or_None)
    Uses PubChem PUG-REST. First gets CID, then MW property, and then synonyms to find a name.
    """
    try:
        cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/cids/JSON"
        cid_resp = requests.get(cid_url, timeout=8)
        cid_resp.raise_for_status()
        cid_data = cid_resp.json()
        cid = cid_data['IdentifierList']['CID'][0]

        # Get molecular weight
        mw_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight/JSON"
        mw_resp = requests.get(mw_url, timeout=8)
        mw_resp.raise_for_status()
        mw_data = mw_resp.json()
        mw = mw_data['PropertyTable']['Properties'][0]['MolecularWeight']

        # Try to fetch synonyms and use the first (often the common name/title)
        syn_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
        syn_resp = requests.get(syn_url, timeout=8)
        preferred_name = None
        if syn_resp.ok:
            try:
                syn_data = syn_resp.json()
                syns = syn_data.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
                if syns:
                    preferred_name = syns[0]
            except Exception:
                preferred_name = None

        return float(mw), cid, preferred_name
    except Exception:
        return None, None, None

# -------------------------
# Session state initialization
# -------------------------
if "mw" not in st.session_state:
    st.session_state["mw"] = 0.0
if "mw_name" not in st.session_state:
    st.session_state["mw_name"] = ""

# Calibration table stored as list of dict rows in session_state
if "cal_rows" not in st.session_state:
    # default three calibration points (mass-based default unit mM is replaced by mg/mL below as we enforce consistency)
    st.session_state.cal_rows = [
        {"conc": 0.10, "conc_unit": "mg/mL", "meas": 0.0000},
        {"conc": 0.20, "conc_unit": "mg/mL", "meas": 0.0000},
        {"conc": 0.50, "conc_unit": "mg/mL", "meas": 0.0000},
    ]

# -------------------------
# Page header
# -------------------------
st.title("⚗️ Chemical Calculator — Version 2.2.1")
st.write("Patch notes: PubChem name display, aligned units, mL defaults, mass/molar concentration support, enforced consistent unit type for dilutions.")

tab1, tab2 = st.tabs(["Chemical Calculator", "Calibration & Dilution Series"])

# -------------------------
# TAB 1: Chemical Calculator
# -------------------------
with tab1:
    st.header("Chemical Calculator")

    # Molecular Weight lookup (PubChem)
    st.subheader("Molecular Weight Lookup")
    col_lookup_name, col_lookup_btn = st.columns([4,1])
    with col_lookup_name:
        compound_name = st.text_input("Compound name (PubChem)", key="lookup_name")
    with col_lookup_btn:
        if st.button("Lookup MW", key="lookup_btn"):
            if compound_name.strip() == "":
                st.warning("Enter a compound name (e.g., glucose, NaCl).")
            else:
                mw_res, cid, pref_name = fetch_molecular_weight_and_name(compound_name.strip())
                if mw_res:
                    st.success(f"Molecular weight: {mw_res:.4f} g/mol")
                    # preferred name display
                    if pref_name:
                        st.info(f"Name found: **{pref_name}**")
                        st.session_state["mw_name"] = pref_name
                    else:
                        st.session_state["mw_name"] = ""
                    st.session_state["mw"] = mw_res
                    st.markdown(f"[View on PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/{cid})")
                else:
                    st.error("Could not find compound on PubChem.")

    st.markdown("---")

    # Mass calculator area
    st.subheader("Mass / Concentration Calculations")

    # First row: MW input and label on right (format 2 decimals)
    c_mw, c_mw_unit = st.columns([4,1])
    with c_mw:
        mw_val = st.number_input("Molecular weight", key="mw_input", value=float(st.session_state.get("mw", 0.0)), format="%.2f")
    with c_mw_unit:
        st.markdown("**g/mol**")

    # Volume + unit compact (default mL)
    c_vol, c_vol_unit = st.columns([3,1])
    with c_vol:
        volume = st.number_input("Volume", key="volume_input", value=1.00, format="%.2f", help="Volume of final solution")
    with c_vol_unit:
        vol_unit = st.selectbox("", ["mL", "L", "µL"], index=0, key="volume_unit_input", label_visibility="collapsed")

    # Concentration + unit compact (allow mass or molar - show two decimals initially)
    c_conc, c_conc_unit = st.columns([3,1])
    with c_conc:
        concentration = st.number_input("Concentration", key="conc_input", value=0.00, format="%.2f", help="Target concentration (choose unit right)")
    with c_conc_unit:
        concentration_unit = st.selectbox("", ["mg/mL", "µg/mL", "ng/mL", "M", "mM", "µM", "nM"], index=0, key="conc_unit_input", label_visibility="collapsed")

    # Calculate required mass button
    if st.button("Calculate Required Mass", key="calc_mass_btn"):
        try:
            # convert final volume to L
            vol_L = float(volume) * (1e-3 if vol_unit == "mL" else (1.0 if vol_unit == "L" else 1e-6))
            # convert concentration to g/L
            if is_molar_unit(concentration_unit) and (mw_val is None or mw_val <= 0):
                st.error("Molecular weight required for molar concentrations.")
            else:
                conc_g_per_L = conc_to_g_per_L(float(concentration), concentration_unit, mw=float(mw_val) if is_molar_unit(concentration_unit) else None)
                # conc_g_per_L is g per L; mass required = conc_g_per_L * vol_L
                mass_g = conc_g_per_L * vol_L
                st.success(f"Mass required: {format_mass(mass_g)}")
                st.caption(f"Breakdown: {concentration} {concentration_unit} × {volume} {vol_unit} -> {mass_g:.6g} g total")
        except Exception as e:
            st.error(f"Calculation error: {e}")

    # Convert to g/L button (shows in g/L)
    st.markdown("Convert concentration → g/L")
    col_gl_left, col_gl_right = st.columns([3,1])
    with col_gl_left:
        pass
    with col_gl_right:
        if st.button("Convert to g/L", key="conv_gl_btn"):
            try:
                if is_molar_unit(concentration_unit) and (mw_val is None or mw_val <= 0):
                    st.error("Molecular weight required for molar concentrations.")
                else:
                    conc_gL = conc_to_g_per_L(float(concentration), concentration_unit, mw=float(mw_val) if is_molar_unit(concentration_unit) else None)
                    st.info(f"Equivalent concentration: {conc_gL:.6f} g/L")
            except Exception as e:
                st.error(f"Error: {e}")

    st.markdown("---")

    # Reverse calculation: mass + volume -> molarity (volume default mL)
    st.subheader("Calculate Molarity from Mass & Volume")

    r_mass_col, r_mass_unit = st.columns([3,1])
    with r_mass_col:
        mass_input = st.number_input("Mass", key="r_mass_input", value=0.00, format="%.2f")
    with r_mass_unit:
        mass_unit = st.selectbox("", ["g", "mg", "µg", "ng"], index=1, key="r_mass_unit_input", label_visibility="collapsed")

    r_vol_col, r_vol_unit = st.columns([3,1])
    with r_vol_col:
        rev_volume = st.number_input("Volume", key="r_volume_input", value=1.00, format="%.2f")
    with r_vol_unit:
        rev_vol_unit = st.selectbox("", ["mL", "L", "µL"], index=0, key="r_volume_unit_input", label_visibility="collapsed")

    if st.button("Calculate Molarity", key="calc_molarity_btn"):
        try:
            # convert mass to grams
            mass_factors = {"g": 1.0, "mg": 1e-3, "µg": 1e-6, "ng": 1e-9}
            m_g = float(mass_input) * mass_factors[mass_unit]
            V_L = float(rev_volume) * (1e-3 if rev_vol_unit == "mL" else (1.0 if rev_vol_unit == "L" else 1e-6))
            mw_float = float(mw_val)
            if mw_float <= 0 or V_L <= 0:
                st.error("Please enter positive MW and volume.")
            else:
                molarity = m_g / mw_float / V_L
                st.success(f"Molarity: {molarity:.6f} mol/L")
        except Exception as e:
            st.error(f"Error: {e}")

# -------------------------
# TAB 2: Calibration & Dilution Series (Option C, Option1 add/remove)
# -------------------------
with tab2:
    st.header("Calibration & Dilution Series — Version 2.2.1")

    # Stock solution definition
    st.subheader("1) Stock solution")
    sc_col_left, sc_col_right = st.columns([3,1])
    with sc_col_left:
        stock_conc = st.number_input("Stock concentration", key="stock_conc_input", value=2.00, format="%.2f", help="Enter concentration of your stock (mass or molar).")
    with sc_col_right:
        stock_conc_unit = st.selectbox("", ["mg/mL", "µg/mL", "ng/mL", "M", "mM", "µM", "nM"], index=0, key="stock_conc_unit_input", label_visibility="collapsed")

    st.caption("Specify how much stock you *have available* (for all calibration tubes combined).")
    stock_vol_col, stock_vol_unit = st.columns([3,1])
    with stock_vol_col:
        stock_avail_volume = st.number_input("Available stock volume", key="stock_avail_input", value=10.00, format="%.2f", help="Total stock volume you have on hand")
    with stock_vol_unit:
        stock_avail_unit = st.selectbox("", ["mL", "L", "µL"], index=0, key="stock_avail_unit_input", label_visibility="collapsed")

    # Final volume per calibration tube
    fv_col, fv_unit_col = st.columns([3,1])
    with fv_col:
        final_vol = st.number_input("Final volume per tube", key="final_vol_input", value=1.00, format="%.2f")
    with fv_unit_col:
        final_vol_unit = st.selectbox("", ["mL", "L", "µL"], index=0, key="final_vol_unit_input", label_visibility="collapsed")

    st.markdown("---")
    st.subheader("2) Calibration points (add / remove rows)")

    # Render each calibration row as a compact horizontal group
    st.write("Enter target concentration and (optional) measured value for each calibration point.")
    rows_container = st.container()
    with rows_container:
        # show header
        header_cols = st.columns([1.2, 2.0, 1.0, 1.0])
        header_cols[0].markdown("**Point**")
        header_cols[1].markdown("**Target conc**")
        header_cols[2].markdown("**Unit**")
        header_cols[3].markdown("**Measured value**")

        for i, row in enumerate(st.session_state.cal_rows):
            a, b, c, d = st.columns([1.2, 2.0, 1.0, 1.0])
            a.write(f"Point {i+1}")
            conc_key = f"cal_conc_{i}"
            unit_key = f"cal_unit_{i}"
            meas_key = f"cal_meas_{i}"
            # Place number + unit on same vertical level by placing them in adjacent columns
            new_conc = b.number_input("", key=conc_key, value=float(row.get("conc", 0.0)), format="%.2f")
            new_unit = c.selectbox("", ["mg/mL", "µg/mL", "ng/mL", "M", "mM", "µM", "nM"], index=(["mg/mL","µg/mL","ng/mL","M","mM","µM","nM"].index(row.get("conc_unit","mg/mL"))), key=unit_key, label_visibility="collapsed")
            new_meas = d.number_input("", key=meas_key, value=float(row.get("meas", 0.0)), format="%.4f")

            # Save back to session_state
            st.session_state.cal_rows[i]["conc"] = new_conc
            st.session_state.cal_rows[i]["conc_unit"] = new_unit
            st.session_state.cal_rows[i]["meas"] = new_meas

    # Global Add/Remove buttons (Option 1): placed below table
    add_col, rem_col = st.columns([1,1])
    with add_col:
        if st.button("Add row", key="add_row_btn"):
            # Default new row uses same unit type as stock if possible; otherwise mass mg/mL
            default_unit = "mg/mL"
            if is_molar_unit(stock_conc_unit):
                default_unit = "mM"
            elif is_mass_unit(stock_conc_unit):
                default_unit = "mg/mL"
            st.session_state.cal_rows.append({"conc": 0.00, "conc_unit": default_unit, "meas": 0.0000})
    with rem_col:
        if st.button("Remove last row", key="remove_row_btn"):
            if len(st.session_state.cal_rows) > 0:
                st.session_state.cal_rows.pop()

    st.markdown("---")

    # Compute dilution plan
    if st.button("Compute Dilution Scheme", key="compute_dilution_btn"):
        try:
            # Enforce consistency: all cal rows must be same type as stock (mass vs molar)
            stock_is_mass = is_mass_unit(stock_conc_unit)
            stock_is_molar = is_molar_unit(stock_conc_unit)
            inconsistent_rows = []
            for i, r in enumerate(st.session_state.cal_rows):
                if stock_is_mass and not is_mass_unit(r["conc_unit"]):
                    inconsistent_rows.append((i+1, r["conc_unit"]))
                if stock_is_molar and not is_molar_unit(r["conc_unit"]):
                    inconsistent_rows.append((i+1, r["conc_unit"]))
            if inconsistent_rows:
                st.error("Unit type mismatch detected. When using Rule B (enforce same type), all calibration points must use the same category of unit as the stock (all mass-based OR all molar).")
                st.write("Mismatches (point, unit):", inconsistent_rows)
                raise ValueError("Unit type mismatch")

            # Convert stock conc to g/L
            mw_float = float(st.session_state.get("mw", 0.0)) if st.session_state.get("mw", 0.0) else None
            if stock_is_molar and (mw_float is None or mw_float <= 0):
                st.error("Molecular weight is required for molar stock concentrations. Please lookup MW in Tab 1.")
                raise ValueError("MW missing for molar stock")

            stock_g_per_L = conc_to_g_per_L(float(stock_conc), stock_conc_unit, mw=mw_float if stock_is_molar else None)
            # Convert volumes to liters
            final_L = float(final_vol) * (1e-3 if final_vol and final_vol_unit == "mL" else (1.0 if final_vol_unit == "L" else 1e-6))
            stock_avail_L = float(stock_avail_volume) * (1e-3 if stock_avail_unit == "mL" else (1.0 if stock_avail_unit == "L" else 1e-6))

            results = []
            total_stock_needed_L = 0.0

            for i, row in enumerate(st.session_state.cal_rows):
                tgt_value = float(row["conc"])
                tgt_unit = row["conc_unit"]
                # Convert target concentration into g/L (must be same type as stock, enforced)
                tgt_g_per_L = conc_to_g_per_L(tgt_value, tgt_unit, mw=mw_float if is_molar_unit(tgt_unit) else None)

                # Use C1*V1 = C2*V2, where C1 = stock_g_per_L, V1 = ? (L)
                if stock_g_per_L <= 0:
                    raise ValueError("Stock concentration must be > 0.")
                V_stock_L = (tgt_g_per_L * final_L) / stock_g_per_L
                V_diluent_L = final_L - V_stock_L
                total_stock_needed_L += V_stock_L

                results.append({
                    "Point": i+1,
                    "Target conc": f"{tgt_value} {tgt_unit}",
                    "Stock volume (mL)": 0.0 if V_stock_L <= 1e-12 else V_stock_L * 1000.0,
                    "Diluent volume (mL)": 0.0 if V_diluent_L <= 1e-12 else V_diluent_L * 1000.0
                })

            res_df = pd.DataFrame(results)
            st.success("Dilution scheme computed")
            st.dataframe(res_df, use_container_width=True)
            st.caption(f"Total stock needed: {total_stock_needed_L*1000:.3f} mL  — Available: {stock_avail_L*1000:.3f} mL")
            if total_stock_needed_L > stock_avail_L + 1e-12:
                st.warning("Warning: Required stock volume exceeds available stock. Increase stock, reduce number of tubes, or prepare a more concentrated stock.")
        except Exception as e:
            if not isinstance(e, ValueError) or "Unit type mismatch" not in str(e):
                st.error(f"Error computing dilution scheme: {e}")

    st.markdown("---")

    # Calibration plot area
    st.subheader("3) Calibration plot (use measured values entered above)")

    plot_col_left, plot_col_right = st.columns([3,1])
    with plot_col_right:
        fit_type = st.selectbox("Fit", ["None", "Linear", "Quadratic"], key="fit_type_input")
        show_equation = st.checkbox("Show fit equation", value=True, key="show_fit_eq")

    # Prepare data for plotting
    concs = []
    meas = []
    for row in st.session_state.cal_rows:
        try:
            val = float(row["conc"])
            if not isfinite(val):
                continue
            concs.append(val)
            meas.append(float(row["meas"]))
        except Exception:
            continue

    if len(concs) >= 1:
        fig, ax = plt.subplots()
        ax.scatter(concs, meas, label="data", zorder=3)

        if fit_type != "None" and len(concs) >= (2 if fit_type == "Linear" else 3):
            deg = 1 if fit_type == "Linear" else 2
            coef = np.polyfit(concs, meas, deg)
            poly = np.poly1d(coef)
            xfit = np.linspace(min(concs), max(concs), 200)
            yfit = poly(xfit)
            ax.plot(xfit, yfit, linestyle="--", label=f"{fit_type} fit")
            if show_equation:
                ymean = np.mean(meas)
                ss_tot = np.sum((meas - ymean) ** 2)
                ss_res = np.sum((np.array(meas) - poly(np.array(concs))) ** 2)
                r2 = 1 - ss_res / ss_tot if ss_tot != 0 else float("nan")
                st.markdown(f"**Fit coefficients:** `{np.array2string(coef, precision=6, separator=', ')}`  •  **R²:** {r2:.4f}")

        ax.set_xlabel("Concentration (entered units)")
        ax.set_ylabel("Measured value (absorbance / signal)")
        ax.grid(True, which="both", linestyle=":", linewidth=0.5)
        ax.legend()
        st.pyplot(fig)
    else:
        st.info("Add at least one calibration point and (optionally) measured values to plot.")

    st.markdown("---")
    st.caption("Tip: All calibration points must use the same unit category as the stock (mass-based OR molar). Use 'Add row' to add points, set final volume (same for each tube), then 'Compute Dilution Scheme'.")
