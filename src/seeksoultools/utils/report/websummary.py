import sys
import json
import base64

class websummary:
    data = {"id": 0}
    def __init__(self, logo, description_json):
        self.data["logo"] = self.encode_image(logo)
        with open(description_json, "r") as f:
            self.data["description"] = json.load(f)

    def encode_image(self, image_path):
        '''
        将图片编码为base64
        :param image_path: 图片路径
        :return: base64编码后的图片
        '''
        with open(image_path, "rb") as f:
            image_data = f.read()
            image_base64 = base64.b64encode(image_data)
            return image_base64.decode()
        
    def get_description(self, key1, key2):
        '''
        根据key1和key2获取描述信息
        :param key1: 描述信息的key1
        :param key2: 描述信息的key2
        :return: 描述信息
        '''
        if key1 in self.data["description"]:
            if key2 in self.data["description"][key1]:
                return self.data["description"][key1][key2]
            else:
                # sys.stderr.write(f"{key2} not found in description.json {key1} section.\n")
                return ""
        else:
            # sys.stderr.write(f"{key1} not found in description.json.\n")
            return ""
    
    def cal_q30(self, key):
        '''
        计算指定key的q30值
        :param key: 指定的key; barcode_q, umi_q
        :return: q30值
        '''
        totol_base = sum([sum(v) for v in self.summary[key].values()])
        base_30 = sum([sum(v[30:]) for v in self.summary[key].values()])
        # 返回指定key的q30值
        return base_30 / totol_base

    def format_comma(self, i: int)->str:
        '''
        格式化整数，返回字符串
        '''
        return f'{i:,}'

    def format_percent(self, f: float)->str:
        '''格式化百分比'''
        if f > 1:
            sys.stderr.write("f is bigger than 1\n")
        return f'{f:.2%}'

    def format_float(self, f: float, n:int=2)->str:
        '''
        格式化浮点数，默认格式为2位小数
        :param f: float类型
        :param n: int类型，默认为2
        :return: str类型
        '''
        return format(f, f'.{n}f')

    def format_for_web(self, data, comp_type, title=""):
        data_json = {
            "id": self.data["id"],
            "name": comp_type,
            "title": title,
            "data": data,
            "description": {
                k: self.get_description(title, k) for k, v in data.items()
            },
        }
        self.data["id"] += 1
        return data_json