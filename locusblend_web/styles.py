"""Streamlit CSS/style layer for the LocusBlend web app.

The stylesheet is injected with a single st.markdown(..., unsafe_allow_html=True)
call, exactly as it was in app.py. Keep the CSS string byte-for-byte: selector
order, media/container queries, comments, and whitespace all matter.
"""

import streamlit as st


def inject_locusblend_css():
    st.markdown(
        """
        <style>
        :root {
            color-scheme: light;
            --lb-bg: #ffffff;
            --lb-surface: #ffffff;
            --lb-sidebar-bg: #f8fafc;
            --lb-text: #111827;
            --lb-muted: #4b5563;
            --lb-border: #e5e7eb;
            --lb-input-bg: #ffffff;
            --lb-input-border: #d1d5db;
        }

        html, body, .stApp {
            background: var(--lb-bg) !important;
            color: var(--lb-text) !important;
            color-scheme: light !important;
        }

        [data-testid="stAppViewContainer"],
        [data-testid="stHeader"],
        [data-testid="stToolbar"] {
            background: var(--lb-bg) !important;
            color: var(--lb-text) !important;
        }

        [data-testid="stSidebar"],
        [data-testid="stSidebarContent"] {
            background: var(--lb-sidebar-bg) !important;
            color: var(--lb-text) !important;
        }

        .block-container {
            background: var(--lb-bg) !important;
            color: var(--lb-text) !important;
        }

        h1, h2, h3, h4, h5, h6,
        p, li, label,
        [data-testid="stMarkdownContainer"],
        [data-testid="stMarkdownContainer"] p,
        [data-testid="stMarkdownContainer"] li,
        [data-testid="stCaptionContainer"],
        [data-testid="stWidgetLabel"],
        [data-testid="stExpander"],
        [data-testid="stExpander"] details,
        [data-testid="stExpander"] summary {
            color: var(--lb-text) !important;
        }

        [data-testid="stCaptionContainer"],
        .stCaptionContainer {
            color: var(--lb-muted) !important;
        }

        input,
        textarea,
        select,
        [data-baseweb="input"] input,
        [data-baseweb="textarea"] textarea,
        [data-baseweb="select"] div,
        [data-baseweb="select"] span {
            background-color: var(--lb-input-bg) !important;
            color: var(--lb-text) !important;
        }

        [data-baseweb="input"],
        [data-baseweb="textarea"],
        [data-baseweb="select"] {
            background-color: var(--lb-input-bg) !important;
            color: var(--lb-text) !important;
        }

        [data-testid="stDataFrame"],
        [data-testid="stTable"],
        [data-testid="stMetric"],
        [data-testid="stMetricLabel"],
        [data-testid="stMetricValue"] {
            color: var(--lb-text) !important;
        }

        [data-testid="stFileUploader"] {
            color: var(--lb-text) !important;
        }

        [data-testid="stFileUploader"] section {
            background: var(--lb-surface) !important;
            color: var(--lb-text) !important;
            border-color: var(--lb-border) !important;
        }

        #MainMenu {
            visibility: hidden !important;
        }

        [data-testid="stToolbar"] {
            display: none !important;
            visibility: hidden !important;
            height: 0 !important;
        }

        [data-testid="stDecoration"] {
            display: none !important;
            visibility: hidden !important;
        }

        [data-testid="stDeployButton"] {
            display: none !important;
            visibility: hidden !important;
        }

        /* Force Streamlit dialogs and modal content to light mode. */
        [data-testid="stDialog"],
        [data-testid="stDialog"] *,
        div[role="dialog"],
        div[role="dialog"] * {
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        [data-testid="stDialog"] button,
        div[role="dialog"] button {
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
        }

        /* Force expanders to remain readable in browser/Streamlit dark mode. */
        [data-testid="stExpander"],
        [data-testid="stExpander"] details,
        [data-testid="stExpander"] summary,
        [data-testid="stExpander"] summary *,
        [data-testid="stExpander"] div,
        [data-testid="stExpander"] p,
        [data-testid="stExpander"] label,
        [data-testid="stExpander"] span {
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Streamlit/BaseWeb controls, menus and dropdown popovers. */
        [data-baseweb="popover"],
        [data-baseweb="popover"] *,
        [data-baseweb="menu"],
        [data-baseweb="menu"] *,
        [data-baseweb="select"],
        [data-baseweb="select"] *,
        [data-baseweb="input"],
        [data-baseweb="input"] *,
        [data-baseweb="textarea"],
        [data-baseweb="textarea"] * {
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Inputs and buttons should remain visible, not hidden. */
        input,
        textarea,
        select,
        button {
            color-scheme: light !important;
        }

        button {
            color: #111827 !important;
        }

        /* Keep Streamlit alerts readable. */
        [data-testid="stAlert"],
        [data-testid="stAlert"] *,
        .stAlert,
        .stAlert * {
            color: #111827 !important;
        }

        /* Keep export/download controls readable. */
        [data-testid="stDownloadButton"],
        [data-testid="stDownloadButton"] *,
        [data-testid="stButton"],
        [data-testid="stButton"] * {
            color: #111827 !important;
        }

        /* --- Force buttons and number steppers to light mode --- */

        /* Normal Streamlit buttons */
        [data-testid="stButton"] button,
        [data-testid="stDownloadButton"] button,
        [data-testid="stFormSubmitButton"] button {
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        /* Streamlit primary buttons, including Update plot / Prepare export when primary */
        button[data-testid="baseButton-primary"],
        [data-testid="baseButton-primary"],
        [data-testid="stButton"] button[data-testid="baseButton-primary"],
        [data-testid="stFormSubmitButton"] button[data-testid="baseButton-primary"] {
            background-color: #ff4b4b !important;
            color: #ffffff !important;
            border: 1px solid #ff4b4b !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        /* Hover states */
        [data-testid="stButton"] button:hover,
        [data-testid="stDownloadButton"] button:hover,
        [data-testid="stFormSubmitButton"] button:hover {
            background-color: #f9fafb !important;
            color: #111827 !important;
            border-color: #9ca3af !important;
        }

        button[data-testid="baseButton-primary"]:hover,
        [data-testid="baseButton-primary"]:hover,
        [data-testid="stButton"] button[data-testid="baseButton-primary"]:hover,
        [data-testid="stFormSubmitButton"] button[data-testid="baseButton-primary"]:hover {
            background-color: #ff3333 !important;
            color: #ffffff !important;
            border-color: #ff3333 !important;
        }

        /* Disabled buttons should be readable, not black */
        [data-testid="stButton"] button:disabled,
        [data-testid="stDownloadButton"] button:disabled,
        [data-testid="stFormSubmitButton"] button:disabled,
        button:disabled,
        button[disabled] {
            background-color: #f3f4f6 !important;
            color: #6b7280 !important;
            border: 1px solid #d1d5db !important;
            opacity: 1 !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        /* BaseWeb number input stepper buttons, including +/- controls */
        [data-baseweb="input"] button,
        [data-baseweb="input"] [role="button"],
        [data-testid="stNumberInput"] button,
        [data-testid="stNumberInput"] [role="button"] {
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #d1d5db !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        /* Icons inside number steppers */
        [data-baseweb="input"] button svg,
        [data-baseweb="input"] [role="button"] svg,
        [data-testid="stNumberInput"] button svg,
        [data-testid="stNumberInput"] [role="button"] svg {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Disabled number stepper buttons */
        [data-baseweb="input"] button:disabled,
        [data-baseweb="input"] [role="button"][aria-disabled="true"],
        [data-testid="stNumberInput"] button:disabled,
        [data-testid="stNumberInput"] [role="button"][aria-disabled="true"] {
            background-color: #f3f4f6 !important;
            color: #9ca3af !important;
            border-color: #d1d5db !important;
            opacity: 1 !important;
        }

        /* Icons inside disabled number steppers */
        [data-baseweb="input"] button:disabled svg,
        [data-baseweb="input"] [role="button"][aria-disabled="true"] svg,
        [data-testid="stNumberInput"] button:disabled svg,
        [data-testid="stNumberInput"] [role="button"][aria-disabled="true"] svg {
            color: #9ca3af !important;
            fill: #9ca3af !important;
            stroke: #9ca3af !important;
        }

        /* Keep selectbox controls light and readable */
        [data-baseweb="select"] > div,
        [data-baseweb="select"] div[role="button"],
        [data-baseweb="select"] svg {
            background-color: #ffffff !important;
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
            color-scheme: light !important;
        }

        /* Checkbox should remain visible on light background */
        [data-testid="stCheckbox"] label,
        [data-testid="stCheckbox"] span,
        [data-testid="stCheckbox"] div {
            color: #111827 !important;
            color-scheme: light !important;
        }

        .lb-sidebar-visual-card {
            container-type: inline-size;
            background: #ffffff !important;
            border: 1px solid var(--lb-border);
            border-radius: 12px;
            padding: 8px;
            margin: 0.25rem 0 0.75rem 0;
            box-shadow: 0 1px 2px rgba(15, 23, 42, 0.06);
            overflow: hidden;
            width: 100%;
            box-sizing: border-box;
        }

        .lb-sidebar-visual-row {
            display: flex;
            flex-direction: row;
            flex-wrap: nowrap;
            align-items: center;
            justify-content: center;
            gap: 8px;
            width: 100%;
            min-width: 0;
            box-sizing: border-box;
        }

        .lb-sidebar-visual-pane {
            min-width: 0;
            max-width: 100%;
            overflow: hidden;
            display: flex;
            align-items: center;
            justify-content: center;
            box-sizing: border-box;
        }

        .lb-sidebar-visual-pane--legend {
            flex: 1.08 1 0;
        }

        .lb-sidebar-visual-pane--diagram {
            flex: 0.92 1 0;
        }

        .lb-sidebar-visual-pane--single {
            flex: 1 1 auto;
        }

        .lb-sidebar-visual-img {
            display: block;
            width: 100%;
            max-width: 100%;
            min-width: 0;
            height: auto;
            max-height: 125px;
            object-fit: contain;
            box-sizing: border-box;
        }

        @container (max-width: 280px) {
            .lb-sidebar-visual-card {
                padding: 6px;
            }

            .lb-sidebar-visual-row {
                gap: 4px;
            }

            .lb-sidebar-visual-img {
                max-height: 115px;
            }
        }

        @container (max-width: 220px) {
            .lb-sidebar-visual-card {
                padding: 4px;
            }

            .lb-sidebar-visual-row {
                gap: 3px;
            }

            .lb-sidebar-visual-img {
                max-height: 100px;
            }
        }

        /* --- Exact Streamlit 1.50+ stBaseButton override --- */

        /* Exact Streamlit button selectors observed in Chrome DevTools.
           Do not target st-emotion-cache-* classes because they are generated. */

        button[data-testid="stBaseButton-secondary"],
        button[data-testid="stBaseButton-tertiary"],
        button[kind="secondary"],
        button[kind="tertiary"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Text, spans, and icons inside secondary/tertiary buttons */
        button[data-testid="stBaseButton-secondary"] *,
        button[data-testid="stBaseButton-tertiary"] *,
        button[kind="secondary"] *,
        button[kind="tertiary"] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Primary Streamlit buttons */
        button[data-testid="stBaseButton-primary"],
        button[kind="primary"] {
            background: #ff4b4b !important;
            background-color: #ff4b4b !important;
            color: #ffffff !important;
            border: 1px solid #ff4b4b !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Text, spans, and icons inside primary buttons */
        button[data-testid="stBaseButton-primary"] *,
        button[kind="primary"] * {
            color: #ffffff !important;
            fill: #ffffff !important;
            stroke: #ffffff !important;
        }

        /* Generic stBaseButton fallback */
        button[data-testid^="stBaseButton"] {
            color-scheme: light !important;
            box-shadow: none !important;
        }

        /* Disabled buttons: make them light gray, not black */
        button[data-testid^="stBaseButton"]:disabled,
        button[data-testid^="stBaseButton"][disabled],
        button[data-testid^="stBaseButton"][aria-disabled="true"],
        button[kind]:disabled,
        button[kind][disabled],
        button[kind][aria-disabled="true"] {
            background: #f3f4f6 !important;
            background-color: #f3f4f6 !important;
            color: #6b7280 !important;
            border: 1px solid #d1d5db !important;
            opacity: 1 !important;
            box-shadow: none !important;
            cursor: not-allowed !important;
            color-scheme: light !important;
        }

        /* Disabled button inner text/icons */
        button[data-testid^="stBaseButton"]:disabled *,
        button[data-testid^="stBaseButton"][disabled] *,
        button[data-testid^="stBaseButton"][aria-disabled="true"] *,
        button[kind]:disabled *,
        button[kind][disabled] *,
        button[kind][aria-disabled="true"] * {
            color: #6b7280 !important;
            fill: #6b7280 !important;
            stroke: #6b7280 !important;
        }

        /* Hover states for enabled secondary/tertiary buttons */
        button[data-testid="stBaseButton-secondary"]:not(:disabled):hover,
        button[data-testid="stBaseButton-tertiary"]:not(:disabled):hover,
        button[kind="secondary"]:not(:disabled):hover,
        button[kind="tertiary"]:not(:disabled):hover {
            background: #f9fafb !important;
            background-color: #f9fafb !important;
            color: #111827 !important;
            border-color: #9ca3af !important;
        }

        /* Hover states for enabled primary buttons */
        button[data-testid="stBaseButton-primary"]:not(:disabled):hover,
        button[kind="primary"]:not(:disabled):hover {
            background: #ff3333 !important;
            background-color: #ff3333 !important;
            color: #ffffff !important;
            border-color: #ff3333 !important;
        }

        /* File uploader buttons and inner labels */
        [data-testid="stFileUploader"] button[data-testid^="stBaseButton"],
        [data-testid="stFileUploader"] button[kind],
        [data-testid="stFileUploader"] button {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        [data-testid="stFileUploader"] button[data-testid^="stBaseButton"] *,
        [data-testid="stFileUploader"] button[kind] *,
        [data-testid="stFileUploader"] button * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* File uploader dropzone and text */
        [data-testid="stFileUploader"] section,
        [data-testid="stFileUploaderDropzone"],
        [data-testid="stFileUploadDropzone"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #e5e7eb !important;
            color-scheme: light !important;
        }

        [data-testid="stFileUploader"] section *,
        [data-testid="stFileUploaderDropzone"] *,
        [data-testid="stFileUploadDropzone"] * {
            color: #111827 !important;
        }

        /* Number input text field */
        [data-testid="stNumberInput"] input,
        [data-baseweb="input"] input {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #d1d5db !important;
            color-scheme: light !important;
        }

        /* Number input stepper controls.
           Streamlit/BaseWeb can use buttons, role=button, or aria-label wrappers. */
        [data-testid="stNumberInput"] button,
        [data-testid="stNumberInput"] button[data-testid^="stBaseButton"],
        [data-testid="stNumberInput"] button[kind],
        [data-testid="stNumberInput"] [role="button"],
        [data-testid="stNumberInput"] div[aria-label],
        [data-testid="stNumberInput"] span[aria-label],
        [data-baseweb="input"] button,
        [data-baseweb="input"] button[data-testid^="stBaseButton"],
        [data-baseweb="input"] button[kind],
        [data-baseweb="input"] [role="button"],
        [data-baseweb="input"] div[aria-label],
        [data-baseweb="input"] span[aria-label] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Number input stepper icons */
        [data-testid="stNumberInput"] button *,
        [data-testid="stNumberInput"] button[data-testid^="stBaseButton"] *,
        [data-testid="stNumberInput"] button[kind] *,
        [data-testid="stNumberInput"] [role="button"] *,
        [data-testid="stNumberInput"] div[aria-label] *,
        [data-testid="stNumberInput"] span[aria-label] *,
        [data-baseweb="input"] button *,
        [data-baseweb="input"] button[data-testid^="stBaseButton"] *,
        [data-baseweb="input"] button[kind] *,
        [data-baseweb="input"] [role="button"] *,
        [data-baseweb="input"] div[aria-label] *,
        [data-baseweb="input"] span[aria-label] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Disabled number input steppers */
        [data-testid="stNumberInput"] button:disabled,
        [data-testid="stNumberInput"] button[disabled],
        [data-testid="stNumberInput"] button[aria-disabled="true"],
        [data-testid="stNumberInput"] [role="button"][aria-disabled="true"],
        [data-baseweb="input"] button:disabled,
        [data-baseweb="input"] button[disabled],
        [data-baseweb="input"] button[aria-disabled="true"],
        [data-baseweb="input"] [role="button"][aria-disabled="true"] {
            background: #f3f4f6 !important;
            background-color: #f3f4f6 !important;
            color: #9ca3af !important;
            border-color: #d1d5db !important;
            opacity: 1 !important;
        }

        /* Disabled number input stepper icons */
        [data-testid="stNumberInput"] button:disabled *,
        [data-testid="stNumberInput"] button[disabled] *,
        [data-testid="stNumberInput"] button[aria-disabled="true"] *,
        [data-testid="stNumberInput"] [role="button"][aria-disabled="true"] *,
        [data-baseweb="input"] button:disabled *,
        [data-baseweb="input"] button[disabled] *,
        [data-baseweb="input"] button[aria-disabled="true"] *,
        [data-baseweb="input"] [role="button"][aria-disabled="true"] * {
            color: #9ca3af !important;
            fill: #9ca3af !important;
            stroke: #9ca3af !important;
        }

        /* Download button */
        [data-testid="stDownloadButton"] button[data-testid^="stBaseButton"],
        [data-testid="stDownloadButton"] button[kind],
        [data-testid="stDownloadButton"] button {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        [data-testid="stDownloadButton"] button[data-testid^="stBaseButton"] *,
        [data-testid="stDownloadButton"] button[kind] *,
        [data-testid="stDownloadButton"] button * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Radio buttons and checkbox text should remain readable */
        [data-testid="stRadio"] *,
        [data-testid="stCheckbox"] * {
            color-scheme: light !important;
        }

        [data-testid="stRadio"] label,
        [data-testid="stRadio"] span,
        [data-testid="stCheckbox"] label,
        [data-testid="stCheckbox"] span {
            color: #111827 !important;
        }

        /* --- Export button and checkbox repair --- */

        /* Repair normal/secondary Streamlit buttons, including Prepare export file. */
        button[data-testid="stBaseButton-secondary"],
        button[kind="secondary"],
        [data-testid="stButton"] button[data-testid="stBaseButton-secondary"],
        [data-testid="stButton"] button[kind="secondary"],
        [data-testid="stFormSubmitButton"] button[data-testid="stBaseButton-secondary"],
        [data-testid="stFormSubmitButton"] button[kind="secondary"],
        [data-testid="stDownloadButton"] button[data-testid="stBaseButton-secondary"],
        [data-testid="stDownloadButton"] button[kind="secondary"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Repair text inside secondary buttons.
           Do not force SVG fill/stroke here; broad SVG rules can break checkbox marks. */
        button[data-testid="stBaseButton-secondary"] span,
        button[kind="secondary"] span,
        [data-testid="stButton"] button[data-testid="stBaseButton-secondary"] span,
        [data-testid="stButton"] button[kind="secondary"] span,
        [data-testid="stFormSubmitButton"] button[data-testid="stBaseButton-secondary"] span,
        [data-testid="stFormSubmitButton"] button[kind="secondary"] span,
        [data-testid="stDownloadButton"] button[data-testid="stBaseButton-secondary"] span,
        [data-testid="stDownloadButton"] button[kind="secondary"] span {
            color: #111827 !important;
        }

        /* Keep primary buttons readable if any remain. */
        button[data-testid="stBaseButton-primary"],
        button[kind="primary"],
        [data-testid="stButton"] button[data-testid="stBaseButton-primary"],
        [data-testid="stButton"] button[kind="primary"],
        [data-testid="stFormSubmitButton"] button[data-testid="stBaseButton-primary"],
        [data-testid="stFormSubmitButton"] button[kind="primary"] {
            background: #ff4b4b !important;
            background-color: #ff4b4b !important;
            color: #ffffff !important;
            border: 1px solid #ff4b4b !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        button[data-testid="stBaseButton-primary"] span,
        button[kind="primary"] span,
        [data-testid="stButton"] button[data-testid="stBaseButton-primary"] span,
        [data-testid="stButton"] button[kind="primary"] span,
        [data-testid="stFormSubmitButton"] button[data-testid="stBaseButton-primary"] span,
        [data-testid="stFormSubmitButton"] button[kind="primary"] span {
            color: #ffffff !important;
        }

        /* Disabled buttons should be light gray and readable. */
        button[data-testid^="stBaseButton"]:disabled,
        button[data-testid^="stBaseButton"][disabled],
        button[data-testid^="stBaseButton"][aria-disabled="true"],
        button[kind]:disabled,
        button[kind][disabled],
        button[kind][aria-disabled="true"] {
            background: #f3f4f6 !important;
            background-color: #f3f4f6 !important;
            color: #6b7280 !important;
            border: 1px solid #d1d5db !important;
            opacity: 1 !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        button[data-testid^="stBaseButton"]:disabled span,
        button[data-testid^="stBaseButton"][disabled] span,
        button[data-testid^="stBaseButton"][aria-disabled="true"] span,
        button[kind]:disabled span,
        button[kind][disabled] span,
        button[kind][aria-disabled="true"] span {
            color: #6b7280 !important;
        }

        /* Checkbox repair:
           Only control the label text. Do NOT override checkbox internal spans/SVGs/marks,
           because Streamlit/BaseWeb uses those to draw the checkmark. */
        [data-testid="stCheckbox"] label,
        [data-testid="stCheckbox"] label p,
        [data-testid="stCheckbox"] label span:not([data-baseweb]) {
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Restore checkbox box visibility without forcing the inner checkmark SVG. */
        [data-testid="stCheckbox"] [data-baseweb="checkbox"] {
            color-scheme: light !important;
        }

        /* Checkbox input should remain selectable and visible. */
        [data-testid="stCheckbox"] input[type="checkbox"] {
            accent-color: #ff4b4b !important;
        }

        /* Remove overly broad checkbox icon overrides from winning.
           This intentionally avoids fill/stroke rules on checkbox descendants. */
        [data-testid="stCheckbox"] svg {
            color: revert !important;
            fill: revert !important;
            stroke: revert !important;
        }

        /* File uploader Browse files button remains secondary-style. */
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-secondary"],
        [data-testid="stFileUploader"] button[kind="secondary"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
        }

        [data-testid="stFileUploader"] button[data-testid="stBaseButton-secondary"] span,
        [data-testid="stFileUploader"] button[kind="secondary"] span {
            color: #111827 !important;
        }

        /* --- Markdown code and documentation readability --- */

        /* Inline markdown code: prevent black background / green text in browser dark mode. */
        [data-testid="stMarkdownContainer"] code,
        [data-testid="stMarkdownContainer"] p code,
        [data-testid="stMarkdownContainer"] li code,
        [data-testid="stMarkdownContainer"] td code,
        [data-testid="stMarkdownContainer"] th code,
        [data-testid="stExpander"] code,
        [data-testid="stExpander"] p code,
        [data-testid="stExpander"] li code,
        code {
            background: #f3f4f6 !important;
            background-color: #f3f4f6 !important;
            color: #111827 !important;
            border: 1px solid #e5e7eb !important;
            border-radius: 4px !important;
            padding: 0.08rem 0.28rem !important;
            font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace !important;
            font-size: 0.88em !important;
            white-space: break-spaces !important;
            color-scheme: light !important;
        }

        /* Code blocks should also stay light. */
        [data-testid="stMarkdownContainer"] pre,
        [data-testid="stMarkdownContainer"] pre code,
        [data-testid="stExpander"] pre,
        [data-testid="stExpander"] pre code,
        pre,
        pre code {
            background: #f8fafc !important;
            background-color: #f8fafc !important;
            color: #111827 !important;
            border-color: #e5e7eb !important;
            color-scheme: light !important;
        }

        /* Markdown tables in the User Guide should remain readable. */
        [data-testid="stMarkdownContainer"] table,
        [data-testid="stMarkdownContainer"] thead,
        [data-testid="stMarkdownContainer"] tbody,
        [data-testid="stMarkdownContainer"] tr,
        [data-testid="stMarkdownContainer"] th,
        [data-testid="stMarkdownContainer"] td {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #e5e7eb !important;
            color-scheme: light !important;
        }

        /* --- Compact sidebar sections --- */

        /* Reduce sidebar inner padding and vertical gaps. */
        [data-testid="stSidebar"] [data-testid="stSidebarContent"] {
            padding-top: 1.0rem !important;
            padding-bottom: 1.0rem !important;
        }

        /* Keep sidebar markdown/caption spacing tighter. */
        [data-testid="stSidebar"] [data-testid="stMarkdownContainer"] {
            margin-bottom: 0.25rem !important;
        }

        [data-testid="stSidebar"] [data-testid="stCaptionContainer"],
        [data-testid="stSidebar"] .stCaptionContainer {
            margin-top: 0.15rem !important;
            margin-bottom: 0.45rem !important;
            line-height: 1.35 !important;
        }

        /* Compact sidebar expanders while keeping them easy to select. */
        [data-testid="stSidebar"] [data-testid="stExpander"] {
            margin-top: 0.25rem !important;
            margin-bottom: 0.45rem !important;
        }

        [data-testid="stSidebar"] [data-testid="stExpander"] details {
            border-radius: 8px !important;
        }

        [data-testid="stSidebar"] [data-testid="stExpander"] summary {
            min-height: 2.35rem !important;
            padding-top: 0.45rem !important;
            padding-bottom: 0.45rem !important;
        }

        /* Reduce oversized blank gaps between Streamlit element containers in sidebar. */
        [data-testid="stSidebar"] [data-testid="stVerticalBlock"] {
            gap: 0.35rem !important;
        }

        [data-testid="stSidebar"] [data-testid="stElementContainer"] {
            margin-bottom: 0.20rem !important;
        }

        /* Make section helper text compact. */
        [data-testid="stSidebar"] p {
            line-height: 1.35 !important;
        }

        /* Keep form submit button close to final section. */
        [data-testid="stSidebar"] [data-testid="stFormSubmitButton"] {
            margin-top: 0.35rem !important;
            margin-bottom: 0.35rem !important;
        }

        /* --- Checkbox and uploader control visibility --- */

        /* Make Streamlit/BaseWeb checkbox boxes and ticks visible in forced light mode. */
        [data-testid="stCheckbox"] [data-baseweb="checkbox"] > div,
        [data-testid="stCheckbox"] [data-baseweb="checkbox"] div[role="checkbox"],
        [data-testid="stCheckbox"] div[role="checkbox"] {
            background-color: #ffffff !important;
            border-color: #ff4b4b !important;
            color-scheme: light !important;
        }

        /* Checked checkbox state. */
        [data-testid="stCheckbox"] input[type="checkbox"]:checked + div,
        [data-testid="stCheckbox"] [aria-checked="true"],
        [data-testid="stCheckbox"] div[role="checkbox"][aria-checked="true"] {
            background-color: #ff4b4b !important;
            border-color: #ff4b4b !important;
        }

        /* Checkbox tick/check icon. Avoid broad rules on all checkbox descendants. */
        [data-testid="stCheckbox"] [aria-checked="true"] svg,
        [data-testid="stCheckbox"] div[role="checkbox"][aria-checked="true"] svg {
            color: #ffffff !important;
            fill: #ffffff !important;
            stroke: #ffffff !important;
        }

        /* Checkbox label remains readable. */
        [data-testid="stCheckbox"] label,
        [data-testid="stCheckbox"] label p,
        [data-testid="stCheckbox"] label span {
            color: #111827 !important;
        }

        /* Browser-native checkbox fallback. */
        [data-testid="stCheckbox"] input[type="checkbox"] {
            accent-color: #ff4b4b !important;
        }

        /* Uploaded file row and remove controls should not appear black in browser dark mode. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"],
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] *,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] * {
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Uploaded-file remove button / small icon buttons. */
        [data-testid="stFileUploader"] button,
        [data-testid="stFileUploader"] button[data-testid^="stBaseButton"],
        [data-testid="stFileUploader"] button[kind],
        [data-testid="stFileUploader"] [role="button"] {
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        [data-testid="stFileUploader"] button svg,
        [data-testid="stFileUploader"] button *,
        [data-testid="stFileUploader"] [role="button"] svg,
        [data-testid="stFileUploader"] [role="button"] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* If Streamlit renders a small dark delete/remove pill, force it light. */
        [data-testid="stFileUploader"] [aria-label*="remove" i],
        [data-testid="stFileUploader"] [aria-label*="delete" i],
        [data-testid="stFileUploader"] [title*="remove" i],
        [data-testid="stFileUploader"] [title*="delete" i] {
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #d1d5db !important;
            color-scheme: light !important;
        }

        [data-testid="stFileUploader"] [aria-label*="remove" i] *,
        [data-testid="stFileUploader"] [aria-label*="delete" i] *,
        [data-testid="stFileUploader"] [title*="remove" i] *,
        [data-testid="stFileUploader"] [title*="delete" i] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* --- Uploaded-file remove-control fix --- */

        /* Uploaded file row containers. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"],
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] > div,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] div,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] > div,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] div {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Uploaded file row text and file-size text. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] span,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] p,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] small,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] span,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] p,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] small {
            color: #111827 !important;
        }

        /* Uploaded-file row buttons, including remove/delete controls. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] button,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] button[data-testid^="stBaseButton"],
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] button[kind],
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] [role="button"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] button,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] button[data-testid^="stBaseButton"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] button[kind],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] [role="button"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            border-radius: 8px !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Uploaded-file row button icons. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] button svg,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] button *,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] [role="button"] svg,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] [role="button"] *,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] button svg,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] button *,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] [role="button"] svg,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] [role="button"] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Directly target common uploader remove/delete button labels. */
        [data-testid="stFileUploader"] button[aria-label*="remove" i],
        [data-testid="stFileUploader"] button[aria-label*="delete" i],
        [data-testid="stFileUploader"] button[aria-label*="clear" i],
        [data-testid="stFileUploader"] button[title*="remove" i],
        [data-testid="stFileUploader"] button[title*="delete" i],
        [data-testid="stFileUploader"] button[title*="clear" i],
        [data-testid="stFileUploader"] [role="button"][aria-label*="remove" i],
        [data-testid="stFileUploader"] [role="button"][aria-label*="delete" i],
        [data-testid="stFileUploader"] [role="button"][aria-label*="clear" i],
        [data-testid="stFileUploader"] [role="button"][title*="remove" i],
        [data-testid="stFileUploader"] [role="button"][title*="delete" i],
        [data-testid="stFileUploader"] [role="button"][title*="clear" i] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            border: 1px solid #d1d5db !important;
            border-radius: 8px !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Icons inside direct remove/delete targets. */
        [data-testid="stFileUploader"] button[aria-label*="remove" i] *,
        [data-testid="stFileUploader"] button[aria-label*="delete" i] *,
        [data-testid="stFileUploader"] button[aria-label*="clear" i] *,
        [data-testid="stFileUploader"] button[title*="remove" i] *,
        [data-testid="stFileUploader"] button[title*="delete" i] *,
        [data-testid="stFileUploader"] button[title*="clear" i] *,
        [data-testid="stFileUploader"] [role="button"][aria-label*="remove" i] *,
        [data-testid="stFileUploader"] [role="button"][aria-label*="delete" i] *,
        [data-testid="stFileUploader"] [role="button"][aria-label*="clear" i] *,
        [data-testid="stFileUploader"] [role="button"][title*="remove" i] *,
        [data-testid="stFileUploader"] [role="button"][title*="delete" i] *,
        [data-testid="stFileUploader"] [role="button"][title*="clear" i] * {
            color: #111827 !important;
            fill: #111827 !important;
            stroke: #111827 !important;
        }

        /* Some Streamlit versions render the file-row action as the last child.
        Keep this scoped to the uploaded file row only. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] > div:last-child,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] > div:last-child *,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] > div:last-child,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] > div:last-child * {
            background-color: #ffffff !important;
            color: #111827 !important;
            border-color: #d1d5db !important;
            color-scheme: light !important;
        }

        /* Do not let the file row action inherit dark theme surfaces. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] [data-baseweb],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] [data-baseweb] {
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* If this still fails, inspect the dark control in DevTools and add its stable data-testid or aria-label selector. Avoid st-emotion-cache-* classes. */

        /* --- Uploaded-file internal scrollbar fix --- */

        /* The remaining black vertical pill in uploaded-file rows is likely an
        internal scrollbar thumb, not a button. Keep this scoped to file uploader. */
        [data-testid="stFileUploader"],
        [data-testid="stFileUploader"] *,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] {
            scrollbar-color: #cbd5e1 #ffffff !important;
            scrollbar-width: thin !important;
            color-scheme: light !important;
        }

        /* WebKit / Chrome scrollbar track inside file uploader. */
        [data-testid="stFileUploader"]::-webkit-scrollbar,
        [data-testid="stFileUploader"] *::-webkit-scrollbar,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"]::-webkit-scrollbar,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] *::-webkit-scrollbar,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"]::-webkit-scrollbar,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] *::-webkit-scrollbar {
            width: 8px !important;
            height: 8px !important;
            background: #ffffff !important;
            background-color: #ffffff !important;
        }

        /* WebKit / Chrome scrollbar thumb inside file uploader. */
        [data-testid="stFileUploader"]::-webkit-scrollbar-thumb,
        [data-testid="stFileUploader"] *::-webkit-scrollbar-thumb,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"]::-webkit-scrollbar-thumb,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] *::-webkit-scrollbar-thumb,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"]::-webkit-scrollbar-thumb,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] *::-webkit-scrollbar-thumb {
            background: #cbd5e1 !important;
            background-color: #cbd5e1 !important;
            border: 2px solid #ffffff !important;
            border-radius: 999px !important;
        }

        /* WebKit / Chrome scrollbar corner inside file uploader. */
        [data-testid="stFileUploader"]::-webkit-scrollbar-corner,
        [data-testid="stFileUploader"] *::-webkit-scrollbar-corner,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"]::-webkit-scrollbar-corner,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] *::-webkit-scrollbar-corner,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"]::-webkit-scrollbar-corner,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] *::-webkit-scrollbar-corner {
            background: #ffffff !important;
            background-color: #ffffff !important;
        }

        /* Keep the uploaded file row surface light even when the scrollbar is present. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"],
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] {
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #111827 !important;
            color-scheme: light !important;
        }

        /* Do not make the file row action black when the browser creates overlay scrollbars. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] *,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] * {
            color-scheme: light !important;
        }

        /* Optional: make uploader file rows less likely to create tiny internal scrollbars. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] {
            overflow: visible !important;
        }

        /* Fallback for Streamlit file-row wrappers that create a tiny scrollable box. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] > div,
        [data-testid="stFileUploader"] [data-testid="stUploadedFile"] > div {
            scrollbar-color: #cbd5e1 #ffffff !important;
            scrollbar-width: thin !important;
            color-scheme: light !important;
        }

        /* --- Uploaded-file delete button fix --- */

        /* DevTools-confirmed target:
        div[data-testid="stFileUploaderDeleteBtn"]
        > button[data-testid="stBaseButton-minimal"][aria-label^="Remove "] */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] {
            background: transparent !important;
            background-color: transparent !important;
            border: none !important;
            box-shadow: none !important;
            color-scheme: light !important;
            display: flex !important;
            align-items: center !important;
            justify-content: center !important;
        }

        /* Exact remove button. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button[data-testid="stBaseButton-minimal"],
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "],
        [data-testid="stFileUploader"] button[kind="minimal"][aria-label^="Remove "] {
            width: 24px !important;
            height: 24px !important;
            min-width: 24px !important;
            min-height: 24px !important;
            max-width: 24px !important;
            max-height: 24px !important;
            padding: 0 !important;
            margin: 0 !important;
            display: inline-flex !important;
            align-items: center !important;
            justify-content: center !important;
            background: #ffffff !important;
            background-color: #ffffff !important;
            color: #64748b !important;
            border: 1px solid #cbd5e1 !important;
            border-radius: 999px !important;
            box-shadow: none !important;
            opacity: 1 !important;
            color-scheme: light !important;
        }

        /* Hover/focus state: keep light, not black. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:hover,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:focus,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:active,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "]:hover,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "]:focus,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "]:active {
            background: #f8fafc !important;
            background-color: #f8fafc !important;
            color: #334155 !important;
            border: 1px solid #94a3b8 !important;
            box-shadow: none !important;
            outline: none !important;
        }

        /* Remove icon SVG sizing and color. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] svg,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "] svg,
        [data-testid="stFileUploader"] button[kind="minimal"][aria-label^="Remove "] svg {
            width: 14px !important;
            height: 14px !important;
            color: #64748b !important;
            fill: none !important;
            stroke: #64748b !important;
            background: transparent !important;
            background-color: transparent !important;
        }

        /* Remove icon path. Critical: do not give the path a white/black background box. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] svg path,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "] svg path,
        [data-testid="stFileUploader"] button[kind="minimal"][aria-label^="Remove "] svg path {
            color: #64748b !important;
            fill: none !important;
            stroke: #64748b !important;
            background: transparent !important;
            background-color: transparent !important;
        }

        /* Hover/focus icon color. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:hover svg,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:hover svg path,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:focus svg,
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] button:focus svg path,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "]:hover svg,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "]:hover svg path {
            color: #334155 !important;
            fill: none !important;
            stroke: #334155 !important;
        }

        /* Keep uploaded file row layout stable. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderFile"] {
            align-items: center !important;
        }

        /* Avoid older broad rules making the delete icon into a dark block. */
        [data-testid="stFileUploader"] [data-testid="stFileUploaderDeleteBtn"] *,
        [data-testid="stFileUploader"] button[data-testid="stBaseButton-minimal"][aria-label^="Remove "] * {
            box-shadow: none !important;
            text-shadow: none !important;
        }

        /* --- WashU Medicine header --- */

        .lb-washu-header {
            width: 100%;
            background: #A51417 !important;
            background-color: #A51417 !important;
            color: #ffffff !important;
            border-radius: 0;
            margin: -0.5rem 0 1.25rem 0;
            padding: 0;
            box-shadow: none;
            color-scheme: light !important;
        }

        .lb-washu-header-inner {
            min-height: 44px;
            display: flex;
            align-items: center;
            justify-content: flex-start;
            padding: 0.45rem 1.25rem;
        }

        .lb-washu-logo {
            display: block;
            height: 30px;
            max-width: 260px;
            object-fit: contain;
        }

        .lb-washu-logo-fallback {
            font-size: 1.25rem;
            font-weight: 700;
            letter-spacing: 0;
            color: #ffffff !important;
        }

        @media (max-width: 700px) {
            .lb-washu-header {
                margin-top: -0.25rem;
                margin-bottom: 1rem;
            }

            .lb-washu-header-inner {
                min-height: 40px;
                padding: 0.4rem 0.85rem;
            }

            .lb-washu-logo {
                height: 25px;
                max-width: 220px;
            }
        }

        /* --- Full-width fixed WashU Medicine header --- */

        :root {
            --lb-washu-header-height: 54px;
            --lb-washu-red: #A51417;
        }

        /* Fixed global header across sidebar + main content. */
        .lb-washu-global-header {
            position: fixed !important;
            top: 0 !important;
            left: 0 !important;
            right: 0 !important;
            width: 100vw !important;
            height: var(--lb-washu-header-height) !important;
            z-index: 999990 !important;
            background: var(--lb-washu-red) !important;
            background-color: var(--lb-washu-red) !important;
            color: #ffffff !important;
            display: flex !important;
            align-items: center !important;
            justify-content: flex-start !important;
            margin: 0 !important;
            padding: 0 !important;
            border: none !important;
            border-radius: 0 !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        /* Header content width. Keep logo aligned with main app content, but header background spans all. */
        .lb-washu-global-inner {
            width: 100% !important;
            height: var(--lb-washu-header-height) !important;
            display: flex !important;
            align-items: center !important;
            justify-content: flex-start !important;
            padding: 0 1.5rem !important;
            box-sizing: border-box !important;
        }

        /* Logo in fixed header. */
        .lb-washu-header-link {
            display: inline-flex !important;
            align-items: center !important;
            text-decoration: none !important;
            color: inherit !important;
        }

        .lb-washu-header-link:visited,
        .lb-washu-header-link:hover,
        .lb-washu-header-link:active {
            text-decoration: none !important;
            color: inherit !important;
        }

        .lb-washu-header-link img {
            display: block !important;
        }

        .lb-washu-global-header .lb-washu-logo {
            display: block !important;
            height: 32px !important;
            max-width: 280px !important;
            object-fit: contain !important;
        }

        /* Fallback text if local logo is unavailable. */
        .lb-washu-global-header .lb-washu-logo-fallback {
            font-size: 1.25rem !important;
            font-weight: 700 !important;
            color: #ffffff !important;
            letter-spacing: 0.01em !important;
        }

        /* Hide/neutralize the legacy non-fixed header container if still rendered. */
        .lb-washu-header {
            display: none !important;
        }

        /* Push the Streamlit app content below the fixed header. */
        [data-testid="stAppViewContainer"] {
            padding-top: var(--lb-washu-header-height) !important;
        }

        /* Push sidebar content below the fixed header. */
        [data-testid="stSidebar"] {
            padding-top: var(--lb-washu-header-height) !important;
        }

        /* Ensure sidebar background begins below header, while header still spans above it. */
        [data-testid="stSidebar"] [data-testid="stSidebarContent"] {
            padding-top: 1rem !important;
        }

        /* Main block should not add another huge top gap. */
        [data-testid="stMain"] .block-container {
            padding-top: 2rem !important;
        }

        /* Streamlit's own top header can otherwise create a blank strip.
           Keep it transparent and visually minimized without removing app controls. */
        [data-testid="stHeader"] {
            background: transparent !important;
            height: 0 !important;
            min-height: 0 !important;
        }

        /* Keep Streamlit sidebar controls selectable above the banner if present. */
        [data-testid="stSidebarCollapseButton"] {
            z-index: 1000000 !important;
        }

        /* Responsive header. */
        @media (max-width: 700px) {
            :root {
                --lb-washu-header-height: 48px;
            }

            .lb-washu-global-inner {
                padding: 0 1rem !important;
            }

            .lb-washu-global-header .lb-washu-logo {
                height: 27px !important;
                max-width: 230px !important;
            }

            [data-testid="stMain"] .block-container {
                padding-top: 1.5rem !important;
            }
        }

        /* --- Keep sidebar expanded for review/demo --- */

        /* Do not attempt to restyle Streamlit's collapsed/reopen button.
        For this review build, prevent users from collapsing the sidebar.
        This avoids the known issue where the sidebar can be hard to reopen. */

        /* Hide only the native "close/collapse sidebar" control when the sidebar is expanded.
        Keep the selector narrow and scoped to the sidebar header. */
        [data-testid="stSidebar"] [data-testid="stSidebarCollapseButton"],
        [data-testid="stSidebar"] [data-testid="stSidebarHeader"] [data-testid="stSidebarCollapseButton"] {
            display: none !important;
            visibility: hidden !important;
            pointer-events: none !important;
        }

        /* Keep the sidebar itself visible and normal. */
        [data-testid="stSidebar"] {
            visibility: visible !important;
            opacity: 1 !important;
        }

        /* Preserve the app-wide WashU header.
        Do not touch native sidebar reopen controls here. */

        /* --- Remove empty sidebar header gap --- */

        /* In this review/demo build, the sidebar collapse button is intentionally hidden.
        Streamlit's sidebar header container then becomes empty but still occupies
        space above the Legend. Collapse only that empty header area. */
        [data-testid="stSidebar"] [data-testid="stSidebarHeader"] {
            height: 0 !important;
            min-height: 0 !important;
            max-height: 0 !important;
            padding: 0 !important;
            margin: 0 !important;
            overflow: hidden !important;
        }

        /* Keep the actual collapse button hidden. */

        /* Do not change the global WashU header or Streamlit main/header layers here. */

        /* --- Sidebar app title --- */

        /* Keep the native Streamlit sidebar collapse mechanics untouched here.
           This build adds an app-level sidebar title instead of trying to restyle
           native sidebar controls. */

        /* Use the sidebar top area for an intentional app title. */
        .lb-sidebar-app-header {
            margin: 0 0 0.85rem 0 !important;
            padding: 0.85rem 0.9rem !important;
            border-radius: 10px !important;
            background: #ffffff !important;
            background-color: #ffffff !important;
            border: 1px solid #e5e7eb !important;
            color: #111827 !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        .lb-sidebar-app-title {
            font-size: 1.05rem !important;
            font-weight: 800 !important;
            line-height: 1.2 !important;
            color: #111827 !important;
            margin: 0 !important;
            padding: 0 !important;
        }

        .lb-sidebar-app-subtitle {
            font-size: 0.78rem !important;
            font-weight: 500 !important;
            line-height: 1.25 !important;
            color: #6b7280 !important;
            margin-top: 0.25rem !important;
            padding: 0 !important;
        }

        /* Reduce only the visible top spacing inside sidebar user content.
           Do not touch native sidebar buttons. */
        [data-testid="stSidebar"] [data-testid="stSidebarUserContent"] {
            padding-top: 0.75rem !important;
        }

        /* Keep the empty Streamlit sidebar header compact. */

        /* Keep the sidebar collapse button hidden for this review/demo build. */

        /* --- Cleaner sidebar top label --- */

        /* Use a compact label instead of a heavy sidebar title card. */
        .lb-sidebar-app-header {
            display: none !important;
        }

        /* Compact sidebar label with a WashU-red accent. */
        .lb-sidebar-app-label {
            margin: 0.25rem 0 1.1rem 0 !important;
            padding: 0.15rem 0 0.15rem 0.7rem !important;
            border-left: 4px solid #A51417 !important;
            background: transparent !important;
            color: #111827 !important;
            box-shadow: none !important;
            color-scheme: light !important;
        }

        .lb-sidebar-app-label span {
            font-size: 0.95rem !important;
            font-weight: 800 !important;
            line-height: 1.2 !important;
            color: #111827 !important;
            letter-spacing: 0.01em !important;
        }

        /* Reduce only the visible top spacing inside sidebar user content.
           Do not touch native sidebar buttons. */
        [data-testid="stSidebar"] [data-testid="stSidebarUserContent"] {
            padding-top: 0.65rem !important;
        }

        /* Keep the empty Streamlit sidebar header compact. */

        /* Keep the sidebar collapse button hidden for this review/demo build. */

        /* --- Top update button hint --- */

        .lb-sidebar-update-hint {
            margin: 0.25rem 0 0.9rem 0 !important;
            padding: 0 !important;
            font-size: 0.76rem !important;
            line-height: 1.25 !important;
            color: #6b7280 !important;
        }

        /* --- Top-level Visualizer / Documentation switcher --- */

        [data-testid="stRadio"] {
            color-scheme: light !important;
        }

        /* Keep the page switcher visually compact. */
        .lb-doc-note {
            color: #6b7280 !important;
            font-size: 0.9rem !important;
            line-height: 1.45 !important;
        }

        /* Documentation page readability. */
        .lb-doc-section {
            margin-top: 1.25rem !important;
            margin-bottom: 1.25rem !important;
        }

        .lb-doc-section h2,
        .lb-doc-section h3 {
            color: #111827 !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )
