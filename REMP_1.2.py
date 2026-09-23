# -*- coding: utf-8 -*-

# To generate permutations
from itertools import product 

#For GUI
from tkinter import Tk, mainloop, LEFT, TOP , scrolledtext, messagebox
from tkinter.ttk import *
from tkinter import *
import tkinter
import tkinter.font as tkFont
import string




#Utility function to convert tuple to string    
def convertTuple(tup): 
    str =  ''.join(tup) 
    return str

'''
Class is defined for Instruction text formatting 
'''
class RichText(Text):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        # Create fonts with cleaner sizing and better readability
        base_font = tkFont.Font(family="Segoe UI", size=11)
        
        bold_font = tkFont.Font(family="Segoe UI", size=11, weight="bold")
        italic_font = tkFont.Font(family="Segoe UI", size=11, slant="italic")
        h1_font = tkFont.Font(family="Segoe UI", size=15, weight="bold")
        normal_font = tkFont.Font(family="Segoe UI", size=11)
        normal_u_font = tkFont.Font(family="Segoe UI", size=11, underline=True)
        b_u_font = tkFont.Font(family="Segoe UI", size=11, weight="bold", underline=True)

        # Configure tags with proper colors and spacing
        self.tag_configure("bold", font=bold_font, foreground="#1a1a1a")
        self.tag_configure("italic", font=italic_font, foreground="#333333")
        self.tag_configure("noraml", font=normal_font, foreground="#1a1a1a")
        self.tag_configure("h1", font=h1_font, foreground="#0066cc", spacing1=10, spacing3=15)
        self.tag_configure("boldunderline", font=b_u_font, foreground="#0066cc")
        self.tag_configure("normalunderline", font=normal_u_font, foreground="#333333")
        
        em = base_font.measure("m")
        lmargin2 = em + base_font.measure("\u2022 ")
        self.tag_configure("bullet", lmargin1=em, lmargin2=lmargin2, font=normal_font)
# Configure hyperlink tag
        link_font = tkFont.Font(family="Segoe UI", size=11, underline=True)
        self.tag_configure("hyperlink", font=link_font, foreground="blue", underline=True)
        self.tag_bind("hyperlink", "<Enter>", self._show_hand_cursor)
        self.tag_bind("hyperlink", "<Leave>", self._show_arrow_cursor)
        self.tag_bind("hyperlink", "<Button-1>", self._click_link)
        
        em = base_font.measure("m")
        lmargin2 = em + base_font.measure("\u2022 ")
        self.tag_configure("bullet", lmargin1=em, lmargin2=lmargin2, font=normal_font)

    def _show_hand_cursor(self, event):
        self.config(cursor="hand2")

    def _show_arrow_cursor(self, event):
        self.config(cursor="")

    def _click_link(self, event):
        index = self.index(f"@{event.x},{event.y}")
        tags = self.tag_names(index)
        for tag in tags:
            if tag.startswith("link_"):
                url = self._links.get(tag)
                if url:
                    import webbrowser
                    webbrowser.open(url)
                break

    def insert_hyperlink(self, index, text, url):
        """Insert a hyperlink at the given index"""
        if not hasattr(self, '_links'):
            self._links = {}
        link_tag = f"link_{len(self._links)}"
        self._links[link_tag] = url
        self.tag_configure(link_tag, foreground="blue", underline=True)
        self.tag_bind(link_tag, "<Enter>", lambda e: self.config(cursor="hand2"))
        self.tag_bind(link_tag, "<Leave>", lambda e: self.config(cursor=""))
        #self.tag_bind(link_tag, "<Button-1>", lambda e, u=url: self._open_url(u))
        self.insert(index, text, (link_tag, "hyperlink"))

    #def _open_url(self, url):
        #import webbrowser
        #webbrowser.open(url)
    def insert_bullet(self, index, text):
        self.insert(index, f"\u2022 {text}\n", "bullet")

class REMP:
    
    def __init__(self, table1,table2):
        self.table1 = table1
        self.table2 = table2
        self.original_seq = ""  # Store original sequence for comparison
        
                
    def generateDegenerateSequences(self,nu_seq):
        
        # Generating triplets from Nucleotide sequence
        nu_seq = nu_seq.translate(str.maketrans('', '', ' \n\t\r')) #Removing White spaces
        nu_seq = nu_seq.upper()
        
        codons_list = {}        
        no_of_codons = int(len(nu_seq)/3)
        nu_seq_trim = nu_seq[:no_of_codons*3]
        self.original_seq = nu_seq_trim  # Store for later comparison
        print("Input Sequence \n",nu_seq_trim)
        # print("no_of_codons", no_of_codons)
        i = 0
        while(i < no_of_codons*3):
            codons_list['condon_'+str(int(i/3))] = [nu_seq[i:i+3]]
            i += 3
        
        # Identifying substitute triplets for given Nucleotide sequence
        for codon in codons_list:
            for triplet in self.table1:
                if codons_list[codon][0] in self.table1[triplet]:
                    for sub in self.table1[triplet] :
                        if sub not in codons_list[codon]:
                            codons_list[codon].append(sub)                            
                  
        #Generating substitute triplets Sequences after substitute triplets       
        codons_substitut_lists = []
        for codon in codons_list:
            codons_substitut_lists.append(codons_list[codon])
        # print("codons_substitut_list",codons_substitut_lists)
        
        # Window contains 4 triplets, restriction site can obtain from 4 triplets 
        window_size = 4
        no_windows = 0
        if no_of_codons < 4:
            no_windows = 1
        else:
            no_windows = no_of_codons - (window_size-1)
        
        # print("no_windows ",no_windows)
               

        
        '''
        Begin Sliding window technique
        
        Window contains 4 triplets , next iteration from previous window right most triplet removed ,
        next triplet after left most of previous window will be added in current window. Will generate the
        all possible degenerate sequence within  window.
        
        Prefix and Postfix of sequence kept same.
        
        '''
        result= {}
        match_count = 0
        for wind in range(0,no_windows):
            window_conodns_list = codons_substitut_lists[wind:wind+4]
            
            # Generating  degenerate sequences within window
            windo_degenerate_sequences = list(product(*window_conodns_list))
            windo_degenerate_sequences_list = [] 
            for tup in windo_degenerate_sequences:
                windo_degenerate_sequences_list.append(convertTuple(tup))
            
            # Identifying the restriction sites in degenerate sequences, each restriction site can obtain from multiple sequences. 
            for Palindrome in self.table2:
                    for seq in windo_degenerate_sequences_list:
                        ind = seq.find(Palindrome)
                        if ind != -1:
                            match_count += 1
                            prefix_ind =  wind*3
                            posfix_ind = prefix_ind+12
                            # Adding prefix and postfix of given sequencre to degenarate sequnces 
                            tot_seq = nu_seq_trim[:prefix_ind]+seq+nu_seq_trim[posfix_ind:]
                            if Palindrome not in result:
                                result[Palindrome] = [tot_seq]                               
                            else:
                                result[Palindrome].append(tot_seq)

        '''
        End Sliding window technique

        '''
        
        #Identifying minimum no of alphabets changed degenerate sequence for each restriction site
        uniq_min_seq = []
        uniq_enzymes = []
        for hexamer in result:
            if hexamer not in uniq_enzymes:
                uniq_enzymes.append(hexamer)
            min_changes = len(nu_seq_trim)
            min_seq = ""
            for degen_seq in result[hexamer]:
                changes_count = 0
                for ind in range(0,len(degen_seq)):
                    if (degen_seq[ind] != nu_seq_trim[ind]):
                        changes_count += 1
                if changes_count <  min_changes:
                    min_changes = changes_count
                    min_seq = degen_seq
                if changes_count == len(hexamer):
                    ind = len(degen_seq)
            if min_seq not in uniq_min_seq:
                uniq_min_seq.append(min_seq)
        
        # Sorting the all degenerate sequences based on minimum no of alphabets changed . 
        for i in range(len(uniq_min_seq)):
            min_ind = -1
            min_changes = len(uniq_min_seq[i])+1
            for j in range(i,len(uniq_min_seq)):
                changes_count = 0
                for ind in range(0,len(uniq_min_seq[j])):
                    if (uniq_min_seq[j][ind] != nu_seq_trim[ind]):
                        changes_count += 1
                if changes_count < min_changes:
                    min_changes = changes_count
                    min_ind = j
            uniq_min_seq[i],uniq_min_seq[min_ind] = uniq_min_seq[min_ind],uniq_min_seq[i]
        
        uniq_seq_enzymes = {}
        
        # Final output for the given input sequence
        for uniq_seq in uniq_min_seq:
            for uniq_enz in uniq_enzymes:
                if uniq_seq.find(uniq_enz) != -1:
                    if uniq_seq not in uniq_seq_enzymes:
                        uniq_seq_enzymes[uniq_seq] = {'enzymes': [(uniq_enz, ' '.join(x+' ('+uniq_enz+')' for x in self.table2[uniq_enz]))], 'original': self.original_seq}
                    else:
                        uniq_seq_enzymes[uniq_seq]['enzymes'].append((uniq_enz, ' '.join(x+' ('+uniq_enz+')' for x in self.table2[uniq_enz])))
        print("Final Output Sequences with Enzymes\n",uniq_seq_enzymes)
        
        return uniq_seq_enzymes
              
       
if __name__=='__main__':

    #Loading Degenerate triplets from file table1.txt 
    table1 = {}    
    with open('./table1.txt', 'r') as document:
        for line in document:
            if line.strip():  
                key, value = line.split(None, 1)
                value = value.strip()
                table1[key] = value.split(',')
    
    #Loading Restriction site and its enzyme from file table2.txt
    table2 = {}
    with open('./table2.txt', 'r') as document:        
        for line in document:
            if line.strip():  
                key, value = line.split(None, 1)  
                value = value.strip()
                table2[key] = value.split(',')
                    
    # print("Table 1 \n",table1)
    # print("Table 2 \n",table2)
          
    #UI Part
    '''
    Instructions regarding  input and output of the program.
    '''
    def clickInstructions():
        
        instru_window = Toplevel(root)
        instru_window.title("Instructions")
        instru_window.geometry("1100x650") 
        bg_color = "ivory3"
        instru_window.configure(background=bg_color)
        
        Tops = Frame(instru_window, width = 800, relief = FLAT,bg = bg_color)
        Tops.pack(side = TOP,expand=True, fill=BOTH)
                
        text = RichText(Tops, width=50, height=15,padx = 70)
        yscrollbar = Scrollbar(Tops, orient=VERTICAL, command=text.yview)
        text['yscroll'] = yscrollbar.set
        yscrollbar.pack(side="right", fill="y")
        text.pack(fill="both", expand=True)
        
        text.insert("end", "                                                   ","noraml")
        text.insert("end", "Instructions", "h1")
        text.insert("end", "                             \n\n","noraml")
        
        text.insert("end", "Note:", "bold")
        text.insert("end", "The files ", "noraml")
        text.insert("end", "table1.txt", "bold")
        text.insert("end", "(Degenerate Codons) and ", "noraml")
        text.insert("end", "table2.txt", "bold")
        text.insert("end", "(Restriction enzyme sites) can be modified by the user.\n\n", "noraml")
        text.insert("end", "Steps involved in using REMP software are explained below using an example:\n", "bold")
        text.insert("end", "To mutate amino acid Glutamate (E) in the following ORF to a Lysine (K).\n","noraml")
        text.insert("end", "GCC GAC TCA GAC CCA TTA ","noraml")
        text.insert("end", "GAA ","bold")
        text.insert("end", "AGC GCC GAC GTG TCA GAC.\n","noraml")
        text.insert("end", " A   D   S   D   P   L  ","noraml")
        text.insert("end", " E   ","bold")
        text.insert("end", "S   A   D   V   S   D  \n\n","noraml")
        text.insert("end","STEP 1\n\n","boldunderline")
        text.insert("end"," Generate a mutant sequence by replacing GAA codon of E with ", "noraml")
        text.insert("end","any codon ", "bold")
        text.insert("end", " of Lysine.\n\n","noraml")
        text.insert("end", "GCC GAC TCA GAC CCA TTA ","noraml")
        text.insert("end", "AAA ","bold")
        text.insert("end", "AGC GCC GAC GTG TCA GAC [OR]\n","noraml")
        text.insert("end", "GCC GAC TCA GAC CCA TTA ","noraml")
        text.insert("end", "AAG ","bold")
        text.insert("end", "AGC GCC GAC GTG TCA GAC \n","noraml")
        text.insert("end", " A   D   S   D   P   L  ","noraml")
        text.insert("end", " K   ","bold")
        text.insert("end", "S   A   D   V   S   D  \n","noraml")
        text.insert("end","____________________________________________________________________________________________________________________\n\n","normal")
        
        text.insert("end","STEP 2\n\n","boldunderline")
        text.insert("end", " Input the mutant sequence to REMP software. In the input sequence, include ","noraml")
        text.insert("end", "complete \n","bold")
        text.insert("end", "codon ","bold")
        text.insert("end", "of the amino acid at 5` end and codons of ","noraml")
        text.insert("end", "at least 4 amino acids on either side","normalunderline")
        text.insert("end", "of intended \n","noraml")
        text.insert("end", "mutation (for efficient amplification in PCR).\n\n","noraml")
        
        text.insert("end", "TCA GAC CCA TTA ","noraml")
        text.insert("end", "AAG ","bold")
        text.insert("end", "AGC GCC GAC GTG [4 codons on sides]\n","noraml")
        text.insert("end", "GAC TCA GAC CCA TTA ","noraml")
        text.insert("end", "AAG ","bold")
        text.insert("end", "AGC GCC GAC GTG TCA [5 codons on sides]\n","noraml")
        text.insert("end", "GCC GAC TCA GAC CCA TTA ","noraml")
        text.insert("end", "AAG ","bold")
        text.insert("end", "AGC GCC GAC GTG TCA GAC [6 codons on sides]\n\n","noraml")
        text.insert("end","Wrong ","bold")
        text.insert("end", "input sequences do not have complete codon of amino acid at 5` end. The codon of amino\n","noraml")
        text.insert("end", "acid at 3` end is not required to be a complete one.\n","noraml")   
        text.insert("end","____________________________________________________________________________________________________________________\n\n","normal")
        
        text.insert("end","STEP 3\n\n","boldunderline")
        text.insert("end", " From the output sequences, generate a mutant primer with a desired restriction site for\n","noraml")
        text.insert("end", "screening the mutation. An ","noraml")
        text.insert("end", "ideal mutant primer ","bold")
        text.insert("end", "is one with silent mutations as ","noraml")
        text.insert("end", "close to intended\n","bold")
        text.insert("end", "mutation ","bold")
        text.insert("end","(E to K in above example) as possible, and, with a ","noraml")
        text.insert("end","least number of base changes ","bold")
        text.insert("end","from\n","noraml")
        text.insert("end","the input sequence, if possible.\n\n","noraml")
        text.insert("end","    - ","bold")
        text.insert("end","REMP output displays mutant sequences carrying a restriction site(s) created via \n","noraml")
        text.insert("end","      silent mutation(s).\n","bold")
        text.insert("end","    - ","bold")
        text.insert("end","The output sequences are sorted according to number of changes from input sequence, with\n","noraml")
        text.insert("end","      the least changes being on ","noraml")
        text.insert("end","top of the list ","bold")
        text.insert("end","and the most at the bottom of the list.\n","noraml")
        text.insert("end","    - ","bold")
        text.insert("end","The name of the restriction enzyme along with its recognition site is displayed at the right\n","noraml")
        text.insert("end","      end of the output sequence.\n\n","noraml")
        text.insert("end","TCAGACCC","noraml")
        text.insert("end","T","boldunderline")
        text.insert("end","TTAAAGAGCGCCGACGTG\t\t\t\t\t\t\t","noraml")
        text.insert("end","DraI (TTTAAA)\n","noraml")
        text.insert("end","TCAGACCCATTAAAGAGCGC","noraml")
        text.insert("end","T","boldunderline")
        text.insert("end","GACGTG\t\t\t\t\t\t\t","noraml")
        text.insert("end","AfeI (AGCGCT)\n\n","noraml")
        text.insert("end","Of the two outputs shown above, both DraI and AfeI are generated via 1 base change, however,\n","noraml")
        text.insert("end","DraI is closest to intended mutation of E(GAA) to K(AAG) and would be ideal for a primer design.\n","noraml")
        text.insert("end","____________________________________________________________________________________________________________________\n\n","normal")
        text.insert("end","NOTES:\n\n","boldunderline")
        text.insert("end"," (1) The output sequence can be copied and final primer with required Tm can be designed\n","noraml")
        text.insert("end"," using a software such as ‘Oligo Calc’ [","noraml")
        text.insert_hyperlink("end","http://biotools.nubic.northwestern.edu/OligoCalc.html.","https://oligocalc.eu/")
        text.insert("end","]\n\n","noraml")
        text.insert("end"," (2) A restriction site generating distinguishable digested DNA pattern of mutant plasmid compared\n","noraml")
        text.insert("end","to a wild type plasmid, can be analyzed using a software such as NEBcutter V2.0 [","noraml")
        text.insert_hyperlink("end","http://nc2.neb.com/NEBcutter2/.","http://nc2.neb.com/NEBcutter2/")
        text.insert("end","]\n\n","noraml")
        text.configure(state ='disabled')

    # Functionality based on input from UI buttons.       
    def doREMP():
        # Show confirmation dialog before running
        confirm = messagebox.askyesno(
            "Confirm Submission",
            "Do the first three nucleotides of your input form a COMPLETE codon of your mRNA ?",
            parent=root
        )
        if not confirm:
            return  # User clicked No, go back to main window

        # Getting input sequence from textbox 
        nu_seq = input_text.get().strip()
        outputText.configure(state ='normal')
        outputText.delete(1.0, END)
        rempObj = REMP(table1,table2)
        uniq_seq_enzymes = rempObj.generateDegenerateSequences(nu_seq)
        xscrollbar.config(command=outputText.xview)
        yscrollbar.config(command=outputText.yview)
        
        # Configure tags BEFORE inserting text - NO FONT SPECIFICATION
        outputText.tag_configure("changed", underline=True, foreground="red")
        
        # Color palette for restriction sites
        colors = ["#FFD700", "#00CED1", "#90EE90", "#FF6B6B", "#87CEEB", "#FFB6C1", "#DDA0DD", "#F0E68C"]
        color_map = {}
        color_idx = 0
        
        # Helper function to blend two colors
        def blend_colors(hex_color1, hex_color2):
            """Blend two hex colors by averaging their RGB values"""
            def hex_to_rgb(hex_color):
                hex_color = hex_color.lstrip('#')
                return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
            
            def rgb_to_hex(rgb):
                return '#{:02x}{:02x}{:02x}'.format(int(rgb[0]), int(rgb[1]), int(rgb[2]))
            
            rgb1 = hex_to_rgb(hex_color1)
            rgb2 = hex_to_rgb(hex_color2)
            blended = tuple((rgb1[i] + rgb2[i]) // 2 for i in range(3))
            return rgb_to_hex(blended)
        
        # Create all site color tags upfront
        for i, color in enumerate(colors):
            tag_name = f"site_color_{i}"
            outputText.tag_configure(tag_name, background=color, foreground="black")
        
        # Cache for merged colors
        merged_tag_cache = {}
        
        if (len(uniq_seq_enzymes) == 0):
            outputText.insert(END, "None")
        else:
            current_line = 1
            for seq in uniq_seq_enzymes:
                # Get data
                original_seq = uniq_seq_enzymes[seq]['original']
                enzyme_list = uniq_seq_enzymes[seq]['enzymes']
                
                # Insert sequence
                outputText.insert(END, seq)
                
                # Find changed positions
                changed_positions = []
                for i in range(min(len(seq), len(original_seq))):
                    if seq[i] != original_seq[i]:
                        changed_positions.append(i)
                
                # Build list of restriction sites and positions
                sites_found = []
                for enzyme_info in enzyme_list:
                    site = enzyme_info[0]
                    if site not in color_map:
                        color_map[site] = colors[color_idx % len(colors)]
                        color_idx += 1
                    
                    # Find all occurrences of this site in sequence
                    search_pos = 0
                    while True:
                        pos = seq.find(site, search_pos)
                        if pos == -1:
                            break
                        sites_found.append((pos, pos + len(site), site, color_map[site]))
                        search_pos = pos + 1
                
                # Map each position to list of colors covering it
                position_colors = {}  # position -> list of colors
                for site_start, site_end, site, site_color in sites_found:
                    for pos in range(site_start, site_end):
                        if pos not in position_colors:
                            position_colors[pos] = []
                        if site_color not in position_colors[pos]:
                            position_colors[pos].append(site_color)
                
                # Apply tags position by position, creating merged colors for overlaps
                for pos in sorted(position_colors.keys()):
                    colors_at_pos = position_colors[pos]
                    
                    if len(colors_at_pos) == 1:
                        # Single color - use existing tag
                        color_idx_for_tag = colors.index(colors_at_pos[0])
                        tag_name = f"site_color_{color_idx_for_tag}"
                    else:
                        # Multiple colors - create merged tag
                        colors_key = tuple(sorted(colors_at_pos))
                        if colors_key not in merged_tag_cache:
                            # Blend all colors together
                            blended = colors_at_pos[0]
                            for c in colors_at_pos[1:]:
                                blended = blend_colors(blended, c)
                            
                            tag_name = f"merged_{blended.replace('#', '')}"
                            if tag_name not in outputText.tag_names():
                                outputText.tag_configure(tag_name, background=blended, foreground="black")
                            merged_tag_cache[colors_key] = tag_name
                        else:
                            tag_name = merged_tag_cache[colors_key]
                    
                    # Apply tag to this position
                    start_idx = f"{current_line}.{pos}"
                    end_idx = f"{current_line}.{pos + 1}"
                    outputText.tag_add(tag_name, start_idx, end_idx)
                
                # Apply changed nucleotide tags
                for pos in changed_positions:
                    start_idx = f"{current_line}.{pos}"
                    end_idx = f"{current_line}.{pos + 1}"
                    outputText.tag_add("changed", start_idx, end_idx)
                
                                # Add enzyme info on same line as sequence with proper spacing
                outputText.insert(END, '\t\t\t\t\t')
                
                # Create enzyme display text and apply colors
                for i, enzyme_info in enumerate(enzyme_list):
                    site = enzyme_info[0]
                    enzyme_display = enzyme_info[1]  # e.g., "PvuII (CAGCTG)"
                    
                    # Add spacing between enzymes
                    if i > 0:
                        outputText.insert(END, '  ')
                    
                    # Get the start position BEFORE inserting (use "end-1c" to exclude trailing newline)
                    enzyme_start_before = outputText.index("end-1c")
                    enzyme_start_line = int(enzyme_start_before.split('.')[0])
                    enzyme_start_col = int(enzyme_start_before.split('.')[1])
                    
                    # Insert the enzyme info
                    outputText.insert(END, enzyme_display)
                    
                    # Find where the site is in the enzyme_display string
                    site_pos_in_display = enzyme_display.find(site)
                    
                    if site_pos_in_display != -1:
                        # Calculate the exact position in the text widget
                        site_col_start = enzyme_start_col + site_pos_in_display
                        site_col_end = site_col_start + len(site)
                        
                        site_start_idx = f"{enzyme_start_line}.{site_col_start}"
                        site_end_idx = f"{enzyme_start_line}.{site_col_end}"
                        
                        # Get the color for this site
                        site_color = color_map.get(site, colors[0])
                        if site_color in colors:
                            color_idx_for_tag = colors.index(site_color)
                            tag_name = f"site_color_{color_idx_for_tag}"
                            
                            # Apply the color tag to the restriction site in enzyme info
                            outputText.tag_add(tag_name, site_start_idx, site_end_idx)
                
                outputText.insert(END, '\n\n')
                current_line += 2  # Move to next sequence line
        
        # NOW disable the text widget
        outputText.configure(state ='disabled')
    
    # Clear the input and outputs
    def Clear():        
        input_text.set("")
        outputText.configure(state ='normal')
        outputText.delete(1.0, END)
        outputText.configure(state ='disabled')

        
    def CopyResults():
        # Enable text widget to read content
        outputText.configure(state='normal')
        
        # Get all text from output
        all_text = outputText.get("1.0", END).strip()
        
        # Copy to clipboard
        root.clipboard_clear()
        root.clipboard_append(all_text)
        
        # Restore disabled state
        outputText.configure(state='disabled')
        
        # Show feedback
        original_text = btnCopy.cget("text")
        btnCopy.config(text="COPIED!")
        root.after(1500, lambda: btnCopy.config(text=original_text))    
    
    # Tkinter GUI
    root = Tk()
    root.title("REMP")
    root.geometry("1200x600")
    bg_color = "gray"
    root.configure(background=bg_color)
    

    Tops = Frame(root, width = 1100, relief = FLAT,bg=bg_color)
    Tops.pack(side = TOP) 
              
    lblTitle = Label(Tops, font = ('Arial', 40, 'bold'), 
                    text = "REMP ", fg = "cyan2", bd = 10, bg=bg_color, justify = LEFT)
    lblTitle.grid(row = 0, column = 0)
    
    lblTitle2 = Label(Tops, font = ('Arial', 20, 'bold'), 
                    text = "REstriction site in Mutant Primer", fg = "cyan2", bd = 10, bg=bg_color, justify = LEFT)
    lblTitle2.grid(row = 0, column = 1,columnspan = 3)    
    
    lblInput = Label(Tops, font = ('Arial', 18), bg=bg_color,
                text = "Mutant Primer",fg = "thistle1", bd = 5, justify=RIGHT)
    lblInput.grid(row = 1, column = 0)
    
    lblInput2_txt = '[Input complete codon of amino acid at 5` end]'
    lblInput = Label(Tops, font = ('Arial', 16), bg=bg_color,
                text = lblInput2_txt,fg = "thistle1", bd = 5, justify=LEFT)
    lblInput.grid(row = 1, column = 1)
    
    input_text = StringVar()

    # Force uppercase and remove spaces automatically
    def validate_input(*args):
        text = input_text.get()
        # Remove spaces/tabs and convert to uppercase
        clean_text = text.replace(" ", "").replace("\t", "").upper()
        # Only update if different to avoid infinite loop
        if text != clean_text:
            input_text.set(clean_text)

    # Trace changes to the input variable (fires on typing AND pasting)
    input_text.trace_add("write", validate_input)

    txtInput = Entry(Tops,  font = ("Arial",16),
                         bd = 2,  width = 90,
                         textvariable = input_text)                          
    txtInput.grid(row = 2, column = 0,columnspan = 5,ipady=4)

    # Bind Enter key to trigger submit
    txtInput.bind('<Return>', lambda event: doREMP())
            

    btnInstru = Button(Tops,  padx = 2, pady = 2, bd = 2, fg = "yellow", 
                        font = ('Arial', 16), width = 12, 
                        text = "INSTRUCTIONS",relief = FLAT,bg = bg_color,
                command=clickInstructions,justify = LEFT).grid(row = 1, column = 4)
    
    
    btnSubmit = Button(Tops, padx = 2, pady = 2, bd = 2, 
                  fg = "green2", font = ('Arial', 16), 
                    width = 8, text = "SUBMIT", bg = bg_color,relief = FLAT,justify = RIGHT,
                    command = doREMP).grid(row = 3, column = 0)

    btnClear = Button(Tops, padx = 2, pady = 2, bd = 2, 
              fg = "OliveDrab1", font = ('Arial', 16), 
                width = 8, text = "CLEAR", bg = bg_color,relief = FLAT,justify = RIGHT,
                command = Clear).grid(row = 3, column = 4)
    
    global btnCopy
    btnCopy = Button(Tops, padx=2, pady=2, bd=2,
                     fg="cyan", font=('Arial', 16),
                     width=8, text="COPY", bg=bg_color, relief=FLAT, justify=RIGHT,
                     command=CopyResults)
    btnCopy.grid(row=3, column=3)
    
    lblOutput = Label(Tops, font = ('Arial', 18), bg=bg_color,
                text = "Mutant Primer With a Restriction Site", bd = 5,fg = "thistle1",justify=LEFT)
    lblOutput.grid(row = 4, column = 0,columnspan = 2)

    xscrollbar = Scrollbar(Tops, orient=HORIZONTAL)
    xscrollbar.grid(row=25, column=0, sticky=N+S+E+W,columnspan = 5)
    
    yscrollbar = Scrollbar(Tops)
    yscrollbar.grid(row=5, column=5, sticky=N+S+E+W,columnspan = 5)
    
    outputText = Text(Tops, wrap=NONE,height = 12, width = 90,font = ("Arial",16),bd = 5, 
                      xscrollcommand=xscrollbar.set,yscrollcommand=yscrollbar.set)      
    outputText.grid(row = 5,column = 0,columnspan = 5)
    outputText.configure(state ='disabled')
    
    if __name__ == "__main__":  
        root.mainloop()
# Suriya 09/09/2026    


