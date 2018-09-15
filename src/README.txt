my_anchor.py:
    arg1: ref_h7n9.fa
    arg2: output.sam
    arg3: eventalign_names.tsv
    
    整个工程的代码，给出所有KF278742.1的reads中各个anchor的评分

integrate_eventalign_KF278742.1_dict.py：
    arg1: eventalign_names.tsv

    实际是my_anchor.py中的一个函数,用以整合eventalign_names.tsv中的信息到一个字典中，并将字典保存到一个文件“integrated_eventalign_KF278742.1_dict.txt”中，成功的话，以后my_anchor.py就不用读eventalign_names.tsv，而是从这个字典文件里得到需要的数据。函数尽可能地计算出打分时需要的所有信息,只留下有用的略去无用的,并以字典形式呈现,使后面的工作能够快捷方便地获取这些数据。问题在于太慢了。
    最终跑出来用了11个小时左右。但是使用eval()函数转换输出结果文件为字典时，出现内存错误，应该是这个字符串太大了。

integrate_eventalign_KF278742.1_format.py:
    arg1: eventalign_names.tsv
    与上面的程序类似，但是输出文件采用自己规定的格式，而不是dict的字符串表示方式，且是增量输出，运行时没有占用巨大内存的dict存在，事实证明好像是比上面快一点。
