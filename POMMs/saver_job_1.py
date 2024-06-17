from upload_to_drive import upload_to_drive
import os

path = "../results/" #path where the file is
file_name = os.listdir(path) # #how it will save the file in google drive-should be equal
file_path = ['../results/' + file for file in file_name] 
folder_id = '1FaChGiDPfpXs3ENXZpkT3iD50zMpA978'

for file, path in zip(file_name, file_path):
    file_id = upload_to_drive(file, path, folder_id)
    print(f'Uploaded file with ID: {file_id}')
