import unittest

class TestMethods(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        
        super(TestMethods, self).__init__(*args, **kwargs)

    def hello(self):
        print("Hello.")

if __name__ == '__main__':
    unittest.main()
