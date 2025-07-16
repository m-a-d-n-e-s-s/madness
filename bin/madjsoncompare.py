import json

class madjsoncompare:
    # """takes two json output files and compares individual keys, accessed by a list of descending keys"""
    def __init__(self, file1, file2):
        """Initialize with two json files, which are read and stored in data1 and data2"""
        if not file1.endswith('.json'):
            raise ValueError("File1 must be a JSON file")
        if not file2.endswith('.json'):
            raise ValueError("File2 must be a JSON file")
        if not file1 or not file2:
            raise ValueError("Both file1 and file2 must be provided")
        if not isinstance(file1, str) or not isinstance(file2, str):
            raise TypeError("file1 and file2 must be strings representing file paths")

        # print out file names
        print('comparing files:', file1, file2)
        self.file1=file1
        self.file2=file2
        self.success=True

        with open(file1, 'r') as f:
            self.data1 = json.load(f)

        with open(file2, 'r') as f:
            self.data2 = json.load(f)

#        print("json1",self.data1)
#        print("json2",self.data2)

    def exitcode(self):
        if self.success:
            return 0
        else:
            return 1

    def compare(self, subsequentkeys, tolerance):
        value1=self.data1
        value2=self.data2
        for i in subsequentkeys:
#            print("key ", i, " in ", value1, value2)
            value1=value1[i]
            value2=value2[i]

        # print("type of key/value",type(i),type(value1))
        success=False
        diff='is different'
        if (type(value1)==float):
            diff=abs(value1-value2)
            success=diff<tolerance
        elif (type(value1)==int):
            diff=0
            success=(value1==value2)
        elif (type(value1)==str):
            success=(value1==value2)
            if success:
                diff=''
        else:
            success=(value1==value2)
            if success:
                diff=''

        if success:
            print("key ", subsequentkeys, " agrees to  ", tolerance,": ", value1, value2, diff)
        else:
            print("key ", subsequentkeys, " differs gt ", tolerance,": ", value1, value2, diff)
        self.success=self.success and success
