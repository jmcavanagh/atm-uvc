"""download all gases at room temperature and pressure (1 atm, 25 C)"""
import fetch
import requests
from bs4 import BeautifulSoup
import pandas as pd
from typing import List


def get_gases_wiki() -> List[str]:
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


def get_gases_chemspi() -> List[str]:
    requests.get("https://api.rsc.org/compounds/v1/filter/properties.json?query=meltingPoint%20%3E%200%20AND%20boilingPoint%20%3C%20100%20AND%20phase%20%3D%20%22gas%22&start=0&count=1000&sort=meltingPoint&order=asc&fields=casRegistryNumber")
    return []


if __name__ == "__main__":
    # gas_list = get_gases_wiki()
    gas_list = get_gases_chemspi()
    for cas_num in gas_list:
        fetch.save_sdf_to_file(fetch.sdf_from_cas(cas_num), "gases_sdf/" + cas_num + ".sdf")
