"""get all CAS numbers from high production volume chemical list highproductionvolume.pdf"""
import PyPDF2
import fetch


def is_int(s: str) -> bool:
    """check if a string is an integer"""
    try:
        int(s)
        return True
    except ValueError:
        return False


def get_cas_numbers() -> list[str]:
    """get all CAS numbers from high production volume chemical list highproductionvolume.pdf"""
    pdf_file_obj = open("highproductionvolume.pdf", "rb")
    pdf_reader = PyPDF2.PdfReader(pdf_file_obj)
    cas_numbers = []

    for page_no in range(len(pdf_reader.pages)):
        page_obj = pdf_reader.pages[page_no]
        page_text = page_obj.extract_text()
        lines = page_text.split("\n")
        for line in lines:
            values = line.split(" ")
            if len(values) >= 2 and is_int(values[0]) and is_int(values[1]):
                s = values[1]
                formatted_s = s[:-3] + "-" + s[-3:-1] + "-" + s[-1]
                cas_numbers.append(formatted_s)

    pdf_file_obj.close()

    return cas_numbers


if __name__ == "__main__":
    cas_nums = get_cas_numbers()
    for cas_num in cas_nums:
        fetch.save_sdf_to_file(fetch.sdf_from_cas(cas_num), "highprod_sdf/" + cas_num + ".sdf")
