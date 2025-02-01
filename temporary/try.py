# Try functionality

# print(bool(1 == 2))

# dict_1 = {
#     (1, 2, 3): "dd"
# }
#
# print(dict_1)

# l_check = ["A", "B"]
#
# def star_check(a: str, b: str):
#     print(a + b)
#
#
# star_check(*l_check)


# @property functionality


# class PropTry:
#
#     def __init__(self):
#         self.a = 2
#
#     @property
#     def static_return(self):
#         return [i for i in range(13)]
#
#     def print_property(self):
#         print(self.static_return)
#
#
# pt1 = PropTry()
# pt1.print_property()

# Environment Variables

import os


os.environ["var_333"] = "Test value 3"
print(f"{os.environ["var_333"]}")
print(f"{os.environ}")
print(f"Get Process ID: {os.getpid()}")
