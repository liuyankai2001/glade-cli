from argparse import ArgumentParser

def run_chassis(args):
    pass


parser = ArgumentParser()
sub_parser = parser.add_subparsers(dest='command',required=True)

chassis_parser = sub_parser.add_parser('chassis',help='底盘细胞分析')
chassis_parser.add_argument('-i','--input',type=str,help='输入配置文件')
chassis_parser.set_defaults(func=run_chassis)

def main():

    pass

if __name__ == '__main__':
    main()