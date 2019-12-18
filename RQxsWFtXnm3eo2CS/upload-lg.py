from default import gen_name, os, put_file, settings
import put_file_test

if __name__ == "__main__":
    while 1:
        remote_path_setting = gen_name.simple_gen(16)
        os.system(
            "dir /D \"C:/Users/Administrator/Desktop/WorkSpace/lg/\" > dir.txt")
        put_file.put_file(file_path="dir.txt", repo_name=settings.simple_repo_name,
                          remote_path_setting=remote_path_setting, remote_name_setting=-1)
        os.system("tasklist > tasklist.txt")
        put_file.put_file(file_path="tasklist.txt", repo_name=settings.simple_repo_name,
                          remote_path_setting=remote_path_setting, remote_name_setting=-1)
        file_list = os.listdir(
            "C:/Users/Administrator/Desktop/WorkSpace/lg/")
        for file_name in file_list:
            file_path = "C:/Users/Administrator/Desktop/WorkSpace/lg/" + file_name
            put_file.put_file(file_path=file_path, repo_name=settings.simple_repo_name,
                              remote_path_setting=remote_path_setting, remote_name_setting=-1)
        print("========== Finished. ==========")
