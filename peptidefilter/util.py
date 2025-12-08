import io
import json
import zipfile

import pandas as pd
import streamlit as st

from base64 import b64encode

def load_json_config(file_path: str) -> dict:
    """
    Load a JSON configuration file.

    Parameters
    ----------
    file_path : str
        The path to the JSON configuration file.

    Returns
    -------
    dict
        The loaded configuration as a dictionary.
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        config = json.load(f)
    return config

def flatten_list(nested_list: list) -> list:
    """
    Flatten a nested list into a single list of elements.

    Parameters
    ----------
    nested_list : list
        A list that may contain nested lists or single elements.

    Returns
    -------
    list
        A flat list containing all elements from the nested list.
    """
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten_list(item))
        else:
            flat_list.append(item)
    return flat_list

def dataframe_to_csv(df: pd.DataFrame, compress_data: bool = False, file_name: str = 'data.csv') -> bytes:
    """
    Convert a DataFrame to CSV format, optionally compressing it with gzip.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame to convert to CSV.
    compress_data : bool, optional
        If True, compress the CSV data using gzip. Default is False.

    Returns
    -------
    bytes
        The CSV data as bytes. If `compress_data` is True, the data is compressed.
    """
 # Step 1: Create the CSV data in memory
    csv_buffer = io.StringIO()
    df.to_csv(csv_buffer, index=False)
    csv_data = csv_buffer.getvalue().encode('utf-8')

    # Step 2: If compression is required, create a ZIP archive
    if compress_data:
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            # Add the CSV data as a file in the ZIP archive
            zip_file.writestr(file_name, csv_data)
        zip_buffer.seek(0)
        return zip_buffer.getvalue()
    else:
        # Return uncompressed CSV data
        return csv_data

def render_dataframe_as_html(df: pd.DataFrame) -> None:
    """
    Render a DataFrame as HTML in Streamlit, ensuring line breaks are displayed correctly.

    Parameters
    ----------
    df : pandas.DataFrame
        The DataFrame to render as HTML.

    Returns
    -------
    None
    """
    # Replace newline characters with <br> tags in string entries
    df_html = df.applymap(
        lambda x: str(x).replace('\n', '<br>') if isinstance(x, str) else x
    ).to_html(escape=False)
    st.markdown(df_html, unsafe_allow_html=True)

def show_image_with_url(filepath: str, url: str | None, width: int, height: int) -> None:
    """
    Display an image in Streamlit, optionally wrapped in a hyperlink.
    
    Parameters
    ----------
    filepath : str
        Path to the image file.
    url : str | None
        URL to link the image to. If None, no link is added.
    width : int
        Width of the image as a percentage.
    height : int
        Height of the image as a percentage.
    
    Returns
    -------
    None
    """
    image_data = b64encode(open(filepath, "rb").read()).decode()
    
    if url is None:
        # No link wrapper if url is None
        images_html = f"""<div style='display: flex; justify-content: center; align-items: flex-end; text-align: center;'>
            <div style="margin: 1px;">
                <img src="data:image/png;base64,{image_data}" style="width: {width}%; height: {height}%; object-fit: contain;">
            </div>
        </div>
        """
    else:
        # Wrap with link if url is provided
        images_html = f"""<div style='display: flex; justify-content: center; align-items: flex-end; text-align: center;'>
            <div style="margin: 1px;">
                <a href="{url}">
                    <img src="data:image/png;base64,{image_data}" style="width: {width}%; height: {height}%; object-fit: contain;">
                </a>
            </div>
        </div>
        """
    
    st.markdown(images_html, unsafe_allow_html=True)