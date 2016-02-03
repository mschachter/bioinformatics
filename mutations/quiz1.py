
from suffix_tree import  *

if __name__ == '__main__':

    st = SuffixTree('panamabananas')

    i = st.find_substring('banana')
    print 'i=',i

    i = st.find_substring('ana')
    print 'i=',i


