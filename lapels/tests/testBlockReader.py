'''
Created on Oct 11, 2012

@author: Shunping Huang
'''


import unittest
import StringIO
from lapels import blockreader

class TestBlockReader(unittest.TestCase):
    def setUp(self):
        f = StringIO.StringIO('''##comment1
##comment2
>block1
1
2
>block2
3
##comment3
>block3
>block4
4
''')
        self.br = blockreader.BlockReader(f)
    
    def test_next(self):
        self.assertEqual(self.br.next(), '>block1\n')
        self.assertEqual(self.br.next(), '1\n')
        self.assertEqual(self.br.next(), '2\n')
        self.assertRaises(StopIteration, self.br.next)
        self.assertRaises(StopIteration, self.br.next)
    
    def test_peek(self):
        self.assertEqual(self.br.peek(), '>block1\n')
        self.assertEqual(self.br.peek(), '>block1\n')
        self.br.next()
        self.br.next()
        self.br.next()        
        self.assertRaises(StopIteration, self.br.next)
        self.assertEqual(self.br.peek(), '>block2\n')
    
    def test_nextBlock(self):
        self.assertTrue(self.br.nextBlock())
        self.assertEqual(self.br.peek(), '>block2\n')
        self.assertTrue(self.br.nextBlock())
        self.assertEqual(self.br.next(), '>block3\n')
        self.assertTrue(self.br.nextBlock())
        self.assertFalse(self.br.nextBlock())


if __name__ == '__main__':
    unittest.main()
