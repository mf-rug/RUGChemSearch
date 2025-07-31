from bs4 import BeautifulSoup
import pandas as pd

# Load your HTML file
with open("/Users/max/Documents/Code/RUGChemSearch/html_save_20250730/Search.html", encoding="utf-8") as f:
    soup = BeautifulSoup(f, "lxml")  # or "html5lib" if lxml fails

# Find the data table by id
table = soup.find("table", {"id": "chemList"})

# Extract rows
rows = []
for tr in table.find_all("tr"):
    cells = [td.get_text(strip=True) for td in tr.find_all("td")]
    if cells:  # skip header or empty rows
        rows.append(cells)

# Define column names manually since header is elsewhere
columns = ["#", "Name", "Formula", "Cas-nr", "Location", "Pot-nr", "Quantity", "Owner"]

# Adjust column selection (some td's are hidden like ID, Permission, etc.)
df = pd.DataFrame(rows, columns=[
    "Index", "ID", "Name", "Formula", "Cas-nr", "Permission", "Location",
    "Pot-nr", "Quantity", "Owner", "OwnerRegNumber", "Request", "Transfer"
])

# Keep only visible/useful columns
df = df[["Name", "Formula", "Cas-nr", "Location", "Pot-nr", "Quantity", "Owner"]]

df.to_csv("/Users/max/Documents/Code/RUGChemSearch/html_save_20250730/table.csv", index=False)

print(df.head())
