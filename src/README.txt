my_anchor.py:
    arg1: ref_h7n9.fa
    arg2: output.sam
    arg3: eventalign_names.tsv
    
    整个工程的代码，给出所有KF278742.1的reads中各个anchor的评分

integrate_eventalign_KF278742.1_dict.py：
    arg1: eventalign_names.tsv

    实际是my_anchor.py中的一个函数,用以整合eventalign_names.tsv中的信息到一个字典中。尽可能地计算出打分时需要的信息,只留下有用的略去无用的,并以字典形式呈现,使后面的工作能够快捷方便地获取这些数据。问题在于太慢了……
