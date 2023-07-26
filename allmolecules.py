"""download all gases at room temperature and pressure (1 atm, 25 C)"""
import fetch
import requests
from bs4 import BeautifulSoup
import pandas as pd
from typing import List


def get_gases() -> List[str]:
    """get table from Wikipedia https://en.wikipedia.org/wiki/List_of_gases"""
    wiki_url = "https://en.wikipedia.org/wiki/List_of_gases"
    table_class = "wikitable sortable"
    response = requests.get(wiki_url)
    soup = BeautifulSoup(response.text, "html.parser")
    table = soup.find("table", {"class": table_class})
    df = pd.read_html(str(table))[0]
    # return column "CAS No" as list
    cas_nums = df["CAS No"].tolist()
    # remove NaN values
    cas_nums = [cas_num for cas_num in cas_nums if str(cas_num) != "nan"]
    return cas_nums


if __name__ == "__main__":
    for cas_num in get_gases():
        fetch.save_sdf_to_file(fetch.sdf_from_cas(cas_num), "gases_sdf/" + cas_num + ".sdf")
