#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# TraceWin Lattice Parsing

**IN DEVELOPMENT**

Comments after `;`


# In[ ]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


from tracewin import parsers


# In[3]:


# These are the elements that are assigned for parsing. 
# Note that there are more elements in tracewin.elements
parsers.ele_parser


# In[5]:


# Read some lattice lines
LINES = open('data/PIP_II_PDR_v_2.dat').readlines()


# In[6]:


eles = []
for line in LINES:
    ele = parsers.parse_lattice_line(line)
    if ele.ele_type == 'end':
        print('end found')
        break
    eles.append(ele)
    
# Complete set
set(ele.ele_type for ele in eles)


# In[7]:


# Parsed elements
eles[0:50]


# # as TraceWin
# 
# Call the `.as_tracewin` method on each element to reproduce the original input.

# In[8]:


for ele in eles[0:50]:
    if ele.ele_type == 'marker':
        continue
    #print(ele, ele.ele_type)
    print(ele.as_tracewin)


# # as ImpactX

# In[9]:


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


# In[10]:


def impactx_lines(eles):
    
    idrift = 0
    iquad = 0
    irfcavity = 0
    isolenoid = 0
    ifieldmap = 0
    isteer = 0
    
    lines = []
    names = []
    freq = 0
    
    for ele in eles:
        t = ele.ele_type
        
        if t == 'drift':
            idrift += 1
            name = f'D{idrift}'
            line = f'{name}.type = drift \n{name}.ds = {ele.L/1000} \n'
            
        elif t == 'quad':
            iquad += 1
            name = f'Q{iquad}'
            line = f'{name}.type = quad \n{name}.ds = {ele.L/1000} \n{name}.k = {ele.G} \n'
            
        elif t == 'field_map':
            if ele.geom == 7700:
                irfcavity += 1
                name = f'RF{irfcavity}'
                line = f'{name}.type = rfcavity \n{name}.ds = {ele.L/1000} \n{name}.escale = {ele.ke} \n{name}.freq = {freq*1e6} \n{name}.phase = {ele.theta} \n{name}.mapsteps = {100} \n'
            elif ele.geom == 10:
                isolenoid += 1
                name = f'SOL{isolenoid}'
                line = f'{name}.type = solenoid_softedge \n{name}.ds = {ele.L/1000}  \n{name}.bscale = {ele.kb} \n{name}.mapsteps = {100} \n'
            else:
                ifieldmap += 1
                name = 'BAD'
                line = ''
                continue
            
        #elif t =='thin_steering':
        #    isteer += 1
        #    name = f'K{isteer}'
        #    line = f'{name}: kicker, L = 0'
        
        elif t == 'freq':
            freq = ele.f
            continue
            
        else:
            name = 'BAD'
            line = ''
            continue
               
         
        # Collect
        lines.append(line)
        names.append(name)
    
    latlines = []
    for c in chunks(names, 8):
        line = ' '.join(c)
        latlines.append(line)
    latline = 'lattice.elements = \n'
    latline += ' \n'.join(latlines)
    latline += ' '

    lines.append(latline)
    
    return lines
            
impactx_lines(eles);


# In[11]:


with open('test.impactx', 'w') as f:
    for line in impactx_lines(eles):
        f.write(line + '\n')


# In[12]:


get_ipython().system('head test.impactx')


# # Fieldmaps
# 
# TODO

# In[49]:


fm = []
g  = []
for ele in eles:
    if ele.ele_type == 'field_map':
        fm.append(ele.filename)
        g.append(ele.geom)


# In[50]:


set(fm)


# In[51]:


set(g)


# In[52]:


gkeys = ('static electric field',
'static magnetic field',
'RF electric field',
'RF magnetic field',
'3D aperture map')
def parse_geom(geom):
    return dict(zip(gkeys, list(reversed(f'{geom:05}'))))

parse_geom(7700)    


# In[53]:


parse_geom(10)


# In[54]:


parse_geom(70)


# In[1]:


# Cleanup
get_ipython().system('rm test.impactx')


# In[ ]:




