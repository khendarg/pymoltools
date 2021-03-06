'''
Copyleft 2019 by Kevin Hendargo, all rights reversed
'''
from __future__ import division, print_function

import pymol
import tempfile
import subprocess
import re
import os
import gzip
import random

import Bio.PDB
try: CODE = Bio.PDB.protein_letters_3to1
except AttributeError: CODE = Bio.PDB.to_one_letter_code

class NumberedSequence(object):
	def __init__(self, iterable):
		self.iterable = iterable
		self.validate()

	def validate(self):
		for i, pair in enumerate(self.iterable):
			if not ((len(pair) == 2) and (type(pair[1]) is int)):
				raise ValueError('Element {} is invalid: {}'.format(i, pair))

	def get_seq(self): return ''.join([pair[0] for pair in self])

	def get_contigs(self):
		rawcontigs = []
		lasti = None
		for pair in self:
			if lasti is None: rawcontigs.append([pair])
			elif pair[1] == (lasti + 1): rawconfigs[-1].append(pair)
			else: rawcontigs.append([pair])
			lasti = pair[1]
		return [NumberedSequence(contig) for contig in rawcontigs]

	def __getitem__(self, index):
		if type(index) is int:
			for pair in self:
				if pair[1] == index: return NumberedSequence([self.iterable[index]])
		else:
			start = self.iterable[0][1] if index.start is None else index.start
			stop = self.iterable[-1][1] if index.stop is None else index.stop
			requested = range(start, stop+step, step)
			substring = []
			for pair in self:
				if pair[1] in requested: substring.append(pair)
			return NumberedSequence(substring)
	def __len__(self): return len(self.iterable)

	def find_closest(self, string, threshold=0.8):
		if len(string) > len(self): raise NotImplementedError
		else:
			idents = {}
			maxident = 0
			for i in range(len(self) - len(string) + 1):
				idents[i] = 0
				for pair, c2 in zip(self.iterable[i:i+len(string)], string):
					if pair[0] == c2: idents[i] += 1
				maxident = max(maxident, idents[i])
				if maxident == len(string): return NumberedSequence(self.iterable[i:i+len(string)])
			for i in idents:
				if idents[i] == maxident: return NumberedSequence(self.iterable[i:i+len(string)])

	def get_range(self):
		start, stop = None, None
		for pair in self:
			if start is None: start = pair[1]
			stop = pair[1]
		return start, stop

	def __iter__(self): return iter(self.iterable)
	def __repr__(self):
		start, stop = self.get_range()
		return '<NumberedSequence range=({}, {}) seq="{}>'.format(start, stop, self.get_seq())

	@staticmethod
	def from_pdb(fn, chain):
		sequence = []
		with open(fn) as f:
			for l in f:
				if l.startswith('ATOM') or l.startswith('HETATM'):
					if 'CA' in l:
						if l[21] != chain: continue
						try: resn = CODE[l[17:20]]
						except KeyError: continue
						try: resi = int(l[22:26])
						except ValueError: print(l)
						sequence.append((resn, resi))
		return NumberedSequence(sequence)

def _resolve_selection(selection=None):
	if selection is None:
		if 'sele' in pymol.cmd.get_names('all', True): selection = 'sele'
		else: selection = 'enabled'
	else: pass

	return selection

def _sele2fa(selection=None):
	tf = tempfile.NamedTemporaryFile()
	pymol.cmd.save(tf.name, selection=_resolve_selection(selection), format='fasta')
	tf.flush()
	tf.seek(0)
	return tf.read().decode('utf-8')

def _sele2seqdict(selection=None):
	tf = tempfile.NamedTemporaryFile()
	pymol.cmd.save(tf.name, selection=_resolve_selection(selection), format='pdb')
	tf.flush()
	tf.seek(0)
	seqdict = {}
	for l in tf: 	
		if l.startswith('ATOM') or (l.startswith('HETATM') and 'MSE' in l): 
			if l[13:15] != 'CA': continue
			resn = CODE[l[17:20]]
			resi = int(l[22:26])
			seqdict[resi] = resn
	return seqdict

def _sele2numbseq(selection=None):
	tf = tempfile.NamedTemporaryFile()
	selection = _resolve_selection(selection)
	pymol.cmd.save(tf.name, selection=selection, format='pdb')
	numbseqs = {}
	for chain in pymol.cmd.get_chains(selection):
		numbseqs[chain] = NumberedSequence.from_pdb(tf.name, chain=chain)
	return numbseqs

def pbcopy(selection=None): 
	'''
DESCRIPTION

	pbcopy copies sequence data to the Mac clipboard

USAGE

	pbcopy [selection]

ARGUMENTS

	selection = str: Object or selection to copy {default:enabled or sele}

SEE ALSO

	xclip
	'''
	s = _sele2fa(selection)
	p = subprocess.Popen(['pbcopy'], stdin=subprocess.PIPE)
	p.communicate(input=s.encode('utf-8'))

def xclip(selection=None): 
	'''
DESCRIPTION

	xclip copies sequence data to the X11 clipboard

USAGE

	xclip [selection]

ARGUMENTS

	selection = str: Object or selection to copy {default:enabled or sele}

SEE ALSO

	pbcopy
	'''
	s = _sele2fa(selection)
	p = subprocess.Popen(['pbcopy'], stdin=subprocess.PIPE)
	p.communicate(input=s.encode('utf-8'))

def _in_spans(n, spans):
	for i, span in enumerate(spans):
		if span[0] <= n <= span[1]: return i
	return None

def select_tmss(selection=None):
	selection = _resolve_selection(selection)
	for objname in pymol.cmd.get_names('objects', True):
		subselection = selection + ' and ' + objname
		for chain in pymol.cmd.get_chains(subselection):
			subsubselection = subselection + ' and c. ' + chain	
			seqdict = _sele2seqdict(subsubselection)
			fasta = '>sequence\n'
			seq = ''.join([seqdict[i] for i in sorted(seqdict)])
			fasta += seq
			p = subprocess.Popen(['hmmtop', '-if=--', '-of=--', '-sf=FAS', '-pi=spred', '-is=pseudo'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
			out, err = p.communicate(input=fasta.encode('utf-8'))
			indices = re.findall('(?:[0-9] *)+$', out)[0].split()[1:]
			if not indices: continue
			indices = [int(x) for x in indices]
			spans = []
			for i in range(0, len(indices), 2):
				spans.append((indices[i], indices[i+1]))

			raw2resi = {}
			resis = sorted(seqdict)
			for i in range(len(seq)):
				raw2resi[i] = resis[i]
			for i, span in enumerate(spans):
				pymol.cmd.select('{}_{}_TM{}'.format(objname, chain, i+1), '{} and c. {} and i. {}-{}'.format(objname, chain, span[0], span[1]))
	pymol.cmd.deselect()
	#print(dir(pymol.cmd))

def spectrum1(x, hstart=240, hstop=0, sstart=60, sstop=60, vstart=90, vstop=90):
	hue = hstart + x * (hstop - hstart)
	saturation = sstart + x * (sstop - sstart)
	value = vstart + x * (vstop - vstart)

	chroma = value * saturation / 100
	hueimage = hue/60
	n = chroma * (1 - abs(hueimage % 2 - 1))
	if 0 <= hueimage < 1: rgbtuple = (chroma, n, 0)
	elif 1 <= hueimage < 2: rgbtuple = (n, chroma, 0)
	elif 2 <= hueimage < 3: rgbtuple = (0, chroma, n)
	elif 3 <= hueimage < 4: rgbtuple = (0, n, chroma)
	elif 4 <= hueimage < 5: rgbtuple = (n, 0, chroma)
	elif 5 <= hueimage < 6: rgbtuple = (chroma, 0, n)
	else: rgbtuple = (0, 0, 0)
	m = value - chroma
	rgb = (
		int((rgbtuple[0] + m)*2.555),
		int((rgbtuple[1] + m)*2.555),
		int((rgbtuple[2] + m)*2.555)
	)
	return '0x%02.x%02.x%02.x' % rgb

def vegalike(i, maxi, gamma=0.0, delta=1.0):
	if maxi <= 9:
		m = i % 9
		if m == 0: rgbtuple = (0x1f, 0x77, 0xb4)
		elif m == 1: rgbtuple = (0xff, 0x7f, 0x0e)
		elif m == 2: rgbtuple = (0x2c, 0xa0, 0x2c)
		elif m == 3: rgbtuple = (0xd6, 0x27, 0x28)
		elif m == 4: rgbtuple = (0x94, 0x67, 0xbd)
		elif m == 5: rgbtuple = (0x8c, 0x56, 0x4b)
		elif m == 6: rgbtuple = (0xe3, 0x77, 0xc2)
		#elif m == 7: rgbtuple = (0x7f, 0x7f, 0x7f)
		elif m == 7: rgbtuple = (0xbc, 0xbd, 0x22)
		elif m == 8: rgbtuple = (0x17, 0xbe, 0xcf)
	#elif maxi <= 20:
	else:
		m = i % 18
		if m == 0: rgbtuple = (0x1f, 0x77, 0xb4)
		elif m == 1: rgbtuple = (0xae, 0xc7, 0xe8)
		elif m == 2: rgbtuple = (0xff, 0x7f, 0x0e)
		elif m == 3: rgbtuple = (0xff, 0xbb, 0x78)
		elif m == 4: rgbtuple = (0x2c, 0xa0, 0x2c)
		elif m == 5: rgbtuple = (0x98, 0xdf, 0x8a)
		elif m == 6: rgbtuple = (0xd6, 0x27, 0x28)
		elif m == 7: rgbtuple = (0xff, 0x98, 0x96)
		elif m == 8: rgbtuple = (0x94, 0x67, 0xbd)
		elif m == 9: rgbtuple = (0xc5, 0xb0, 0xd5)
		elif m == 10: rgbtuple = (0x8c, 0x56, 0x4b)
		elif m == 11: rgbtuple = (0xc4, 0x9c, 0x94)
		elif m == 12: rgbtuple = (0xe3, 0x77, 0xc2)
		elif m == 13: rgbtuple = (0xf7, 0xb6, 0xd2)
		#elif m == 14: rgbtuple = (0x7f, 0x7f, 0x7f)
		#elif m == 14: rgbtuple = (0xc7, 0xc7, 0xc7)
		elif m == 14: rgbtuple = (0xbc, 0xbd, 0x22)
		elif m == 15: rgbtuple = (0xdb, 0xdb, 0x8d)
		elif m == 16: rgbtuple = (0x17, 0xbe, 0xcf)
		elif m == 17: rgbtuple = (0x9e, 0xda, 0xe5)

	def _desaturate(rgbtuple, x):
		mean = sum(rgbtuple)/3.
		#if x <= 0: return rgbtuple
		#elif x >= 1: return tuple([int(mean) for i in range(3)])
		#else: return tuple([int(rgbtuple[i]*(1-x) + mean*(x)) for i in range(3)])

		#this should allow for increasing saturation (use negative x)
		return tuple([int(rgbtuple[i]*(1-x) + mean*(x)) for i in range(3)])

	def _mult(rgbtuple, x):
		return tuple([int(rgbtuple[i]*x) for i in range(3)])
		
	rgbtuple = _desaturate(rgbtuple, gamma)
	rgbtuple = _mult(rgbtuple, delta)

	#if gamma != 1.0:
	#	rgbtuple = tuple([int(0xff * (x/0xff)**gamma) for x in rgbtuple])

	return '0x%02.x%02.x%02.x' % rgbtuple

def paint_pfam(selection=None, gray=True):
	gray = _str2bool(gray)
	selection = _resolve_selection(selection)

	dom2clan = {}
	if 'PFAMCLANSDB' not in os.environ: raise KeyError('Environment variable $PFAMCLANSDB unset')
	with gzip.GzipFile(filename=os.environ['PFAMCLANSDB']) as f:
		for l in f:
			if l.strip():
				sl = l.split('\t')
				dom2clan[sl[0]] = sl[1]

	domains = {}
	clans = {}

	for objname in pymol.cmd.get_names('objects', True):
		subselection = selection + ' and ' + objname
		for chain in pymol.cmd.get_chains(subselection):
			subsubselection = subselection + ' and c. ' + chain

			if gray: pymol.cmd.color('gray', subsubselection)

			seqdict = _sele2seqdict(subsubselection)
			numbseq = _sele2numbseq(subsubselection)[chain]
			fasta = '>sequence\n' + numbseq.get_seq()

			if 'PFAMDB' not in os.environ: raise KeyError('Environment variable $PFAMDB unset')

			raw_tf = tempfile.NamedTemporaryFile()
			pymol.cmd.save(raw_tf.name, selection=subsubselection, format='fasta')
			raw_tf.flush()
			raw_tf.seek(0)

			tf = tempfile.NamedTemporaryFile()
			for l in raw_tf: tf.write(l.replace('?', 'X'))
			tf.flush()

			p = subprocess.Popen(['hmmscan', '--cpu', '4', '--noali', '--cut_ga', '-o', '/dev/null', '--domtblout', '/dev/stdout', os.environ['PFAMDB'], tf.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

			out, err = p.communicate()

			for l in out.split('\n'):
				if not l.strip(): continue
				elif l.startswith('#'): continue
				else:
					sl = l.split()
					domshort = sl[0]
					domacc = sl[1]
					domaccshort = sl[1].split('.')[0]
					domexpect = float(sl[6])
					envspan = [int(x) for x in sl[19:21]]
					domsel = subsubselection + ' and i. {}-{}'.format(*envspan)

					domain = {'acc':domacc, 'short':domshort, 'accshort':domaccshort, 'span':envspan, 'selection':domsel, 'expect':domexpect}

					try: domains[domaccshort].append(domain)
					except KeyError: domains[domaccshort] = [domain]

					domclan = dom2clan[domaccshort]
					
					try: clans[domclan].append(domain)
					except KeyError: clans[domclan] = [domain]

		for i, clanacc in enumerate(sorted(clans)):
			try: x = i/(len(clans)-1)
			except ZeroDivisionError: x = 0.5
			for j, domain in enumerate(clans[clanacc]):
				
				saturation = random.randint(5, 8)*10
				value = random.randint(4, 9)*10
				#def spectrum1(x, hstart=240, hstop=0, sstart=60, sstop=60, vstart=90, vstop=90):
				
				color = spectrum1(x=x, sstart=saturation, sstop=saturation, vstart=value, vstop=value)
				color = vegalike(i, len(clans), gamma=0.25 - j/6., delta=1.2*0.8**j)

				pymol.cmd.color(color, domain['selection'])

	for clanname in clans:
		if clanname: 
			print('List of domains in clan {}:'.format(clanname))
			print('===============================')
		else: 
			print('List of domains in clan Lowborn:')
			print('================================')
		clselector = ''
		for domain in clans[clanname]:
			print('{} ({:0.0e}): {} ({})'.format(domain['acc'], domain['expect'], domain['short'], domain['selection']))
			clselector += '({}) '.format(domain['selection'])
			
		if clanname: pymol.cmd.select(clanname, clselector)
		else: pymol.cmd.select('Clanless', clselector)
		for domain in clans[clanname]:
			pymol.cmd.select(domain['short'], domain['selection'])
		print()
		
	#for domname in domains:
	#	selector = ''
	#	for domain in domains[domname]:
	#		selector += '({}) '.format(domain['selection'])

	#		print(

	#	pymol.cmd.select(domname, selector)
	pymol.cmd.deselect()
		



def paint_tmss(selection=None, gray=False):
	gray = _str2bool(gray)
	selection = _resolve_selection(selection)

	for objname in pymol.cmd.get_names('objects', True):
		subselection = selection + ' and ' + objname
		for chain in pymol.cmd.get_chains(subselection):
			subsubselection = subselection + ' and c. ' + chain
			seqdict = _sele2seqdict(subsubselection)
			numbseq = _sele2numbseq(subsubselection)[chain]
			fasta = '>sequence\n' + numbseq.get_seq()

			p = subprocess.Popen(['hmmtop', '-if=--', '-of=--', '-sf=FAS', '-pi=spred', '-is=pseudo'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
			out, err = p.communicate(input=fasta.encode('utf-8'))
			indices = re.findall('(?:[0-9] *)+$', out)[0].split()[1:]
			seq = numbseq.get_seq()
			if not indices: continue
			indices = [int(x) for x in indices]
			spans = []
			paintme = []
			for i in range(0, len(indices), 2):
				span = seq[indices[i]-1:indices[i+1]]
				paintme.append(numbseq.find_closest(span))

			if gray: pymol.cmd.color('gray', subsubselection)
			for i, span in enumerate(paintme):
				try: x = 1. * i / (len(paintme)-1)
				except ZeroDivisionError: x = 0.5

				start, stop = span.get_range()
				pymol.cmd.color(spectrum1(x), subsubselection + ' and i. {}-{}'.format(start, stop))

def chview(view='cartoon', selection=None, hidelig=False):
	selection = _resolve_selection(selection)
	if not hidelig: selection += ' and not het'
	pymol.cmd.hide('everything', selection)
	pymol.cmd.show(view, selection)

def _str2bool(s):
	if not s: return False
	elif type(s) is int: return bool(s)
	elif type(s) is float: return bool(s)
	elif type(s) is str:
		if s.lower().startswith('f'): return False
		elif s.lower().startswith('n'): return False
		elif s.lower().startswith('0'): return False

		elif s.lower().startswith('t'): return True
		elif s.lower().startswith('y'): return True
		elif s.lower().startswith('1'): return True
	else: return bool(s)

def rasmolcolor(selection=None):
	selection = _resolve_selection(selection)

	pymol.cmd.color('0xE60A0A', selection + ' and r. ASP+GLU')
	pymol.cmd.color('0xE6E600', selection + ' and r. CYS+MET')
	pymol.cmd.color('0x145AFF', selection + ' and r. LYS+ARG')
	pymol.cmd.color('0xFA9600', selection + ' and r. SER+THR')
	pymol.cmd.color('0x3232AA', selection + ' and r. PHE+TYR')
	pymol.cmd.color('0x00DCDC', selection + ' and r. ASN+GLN')
	pymol.cmd.color('0xEBEBEB', selection + ' and r. GLY')
	pymol.cmd.color('0x0F820F', selection + ' and r. LEU+VAL+ILE')
	pymol.cmd.color('0xC8C8C8', selection + ' and r. ALA')
	pymol.cmd.color('0xB45AB4', selection + ' and r. TRP')
	pymol.cmd.color('0x8282D2', selection + ' and r. HIS')
	pymol.cmd.color('0xDC9682', selection + ' and r. PRO')

def shapelycolor(selection=None):
	selection = _resolve_selection(selection)

	pymol.cmd.color('0xA00042', selection + ' and r. ASP+THR')
	pymol.cmd.color('0x660000', selection + ' and r. GLU')
	pymol.cmd.color('0xFFFF70', selection + ' and r. CYS')
	pymol.cmd.color('0xB8A042', selection + ' and r. MET+TYR')
	pymol.cmd.color('0x4747B8', selection + ' and r. LYS')
	pymol.cmd.color('0x00007C', selection + ' and r. ARG')
	pymol.cmd.color('0xFF4C4C', selection + ' and r. SER+GLN')
	pymol.cmd.color('0x534C42', selection + ' and r. PHE+PRO+TRP')
	pymol.cmd.color('0xFF7C70', selection + ' and r. ASN')
	pymol.cmd.color('0xFFFFFF', selection + ' and r. GLY+VAL')
	pymol.cmd.color('0x004C00', selection + ' and r. ILE')
	pymol.cmd.color('0x455E45', selection + ' and r. LEU')
	pymol.cmd.color('0x8CFF8C', selection + ' and r. ALA')
	pymol.cmd.color('0x7070FF', selection + ' and r. HIS')

def memecolor(selection=None): 
	selection = _resolve_selection(selection)

	pymol.cmd.color('blue', selection + ' and r. ALA+CYS+PHE+ILE+LEU+MET+VAL+TRP')
	pymol.cmd.color('magenta', selection + ' and r. ASP+GLU')
	pymol.cmd.color('orange', selection + ' and r. GLY')
	pymol.cmd.color('pink', selection + ' and r. HIS')
	pymol.cmd.color('red', selection + ' and r. LYS+ARG')
	pymol.cmd.color('green', selection + ' and r. ASN+GLN+SER+THR')
	pymol.cmd.color('yellow', selection + ' and r. PRO')
	pymol.cmd.color('cyan', selection + ' and r. TYR')

def clustalcolor(selection=None): 
	selection = _resolve_selection(selection)

	pymol.cmd.color('blue', selection + ' and r. ALA+PHE+ILE+LEU+MET+VAL+TRP')
	pymol.cmd.color('pink', selection + ' and r. CYS')
	pymol.cmd.color('magenta', selection + ' and r. ASP+GLU')
	pymol.cmd.color('orange', selection + ' and r. GLY')
	pymol.cmd.color('red', selection + ' and r. LYS+ARG')
	pymol.cmd.color('green', selection + ' and r. ASN+GLN+SER+THR')
	pymol.cmd.color('yellow', selection + ' and r. PRO')
	pymol.cmd.color('teal', selection + ' and r. HIS+TYR')


pymol.cmd.extend('rasmolcolor', rasmolcolor)
pymol.cmd.extend('shapelycolor', shapelycolor)
pymol.cmd.extend('memecolor', memecolor)
pymol.cmd.extend('clustalcolor', clustalcolor)

pymol.cmd.extend('pbcopy', pbcopy)
pymol.cmd.extend('xclip', xclip)

pymol.cmd.extend('cv', chview)

pymol.cmd.extend('tmselect', select_tmss)
pymol.cmd.extend('tmpaint', paint_tmss)

pymol.cmd.extend('dompaint', paint_pfam)

