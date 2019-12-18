import requests
import random
import string
import time
import hashlib


def gen(length: int):
    ret = ""
    for _ in range(length//8):
        ret += ''.join(random.sample(string.digits + string.ascii_letters, 8))
    return ret


if __name__ == "__main__":
    random.seed(time.time_ns())
    while (1):
        random_string = gen(131072)
        data = {"text": random_string}
        result = requests.post(
            url="https://skywt.cn/paste/api/push.php", data=data)
        print(result.text)
