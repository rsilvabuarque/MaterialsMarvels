from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
import subprocess
from time import sleep
import requests

MAX_DELAY = 5
BUTTON_SELECTOR = 'jp-button[data-command="notebook:run-cell-and-select-next"]'
CHROMEDRIVER_PATH = "/Users/sachin/Documents/Coding/Self-Assembly/api/chromedriver" # Executable varies based on type of OS and version of Chrome
TOKEN = "f9a3bd4e9f2c3be01cd629154cfb224c2703181e050254b5" # Arbitrary token chosen for authentication

# Run Jupyter Notebook Server
p = subprocess.Popen(["python3", "-m", "jupyter", "notebook", "--port=8888", f"--IdentityProvider.token={TOKEN}", "--no-browser"])
sleep(5) # Wait for server to start

# Initialize Driver
options = Options()
# options.add_argument("--headless")
driver = webdriver.Chrome(options=options, service=Service(CHROMEDRIVER_PATH))

# Fetch Jupyter Webpage
driver.get(f"http://localhost:8888/notebooks/export_visuals.ipynb?token={TOKEN}")

# Run All Cells
try:
    # Wait for the Run button to load
    button = WebDriverWait(driver, MAX_DELAY).until(EC.presence_of_element_located((By.CSS_SELECTOR, BUTTON_SELECTOR)))

    # Click the Run button to run first cell
    button.click()

    # Wait for visualization to appear
    WebDriverWait(driver, MAX_DELAY).until(EC.presence_of_element_located((By.TAG_NAME, "canvas")))
    sleep(1)

    # Click the button to run second cell
    button.click()
    sleep(1)
except TimeoutException:
    print("Loading took too much time!")

# Kill Jupyter notebook process
p.kill()
# Shutdown Jupyter notebook server
requests.post(f"http://localhost:8888/api/shutdown?token={TOKEN}")