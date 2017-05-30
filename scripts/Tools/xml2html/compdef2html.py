#!/usr/bin/env python

"""Generator of html file for compsets
"""

# Typically ignore this.
# pylint: disable=invalid-name

# Disable these because this is our standard setup
# pylint: disable=wildcard-import,unused-wildcard-import,wrong-import-position

import os, sys, re, glob
import datetime
import collections

CIMEROOT = os.environ.get("CIMEROOT")
if CIMEROOT is None:
    raise SystemExit("ERROR: must set CIMEROOT environment variable")
sys.path.append(os.path.join(CIMEROOT, "scripts", "Tools"))

SRCROOT = os.environ.get("SRCROOT")
if SRCROOT is None:
    SRCROOT = os.path.dirname(CIMEROOT)
os.environ['SRCROOT'] = SRCROOT

from standard_script_setup import *
from CIME.utils            import expect, get_model, get_cime_root
from CIME.XML.entry_id     import GenericXML
from CIME.XML.files        import Files
from CIME.XML.component    import Component
from CIME.XML.compsets     import Compsets
from CIME.case             import Case
import xml.etree.ElementTree as etree

# check for  dependency module
try:
    import jinja2
except:
    raise SystemExit("ERROR: compdef2html.py depends on the jinja2 template module. " /
                     "Install using 'pip --user install jinja2'")

# global variables
_now = datetime.datetime.now().strftime('%Y-%m-%d')
logger = logging.getLogger(__name__)
_cime_comps = {'data_comps' : ['datm','desp','dice','dlnd','docn','drof','dwav'],
               'stub_comps' : ['satm','sesp','sice','slnd','socn','srof','swav'],
               'xcpl_comps' : ['xamp','xexp','xice','xlnd','xocn','xrof','xwav'] }

###############################################################################
def commandline_options():
###############################################################################

    """ Process the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Read all the config_compset.xml files and generate a corresponding HTML file.')

    CIME.utils.setup_standard_logging_options(parser)

    parser.add_argument('--htmlfile', nargs=1, required=True,
                        help='Fully quailfied path to output HTML file.')

    parser.add_argument('--version', nargs=1, required=True,
                        help='Model version (e.g. CESM2.0)')

    options = parser.parse_args()

    CIME.utils.parse_args_and_handle_standard_logging_options(options)

    return options

###############################################################################
def get_desc(compset, desc_file):
###############################################################################

    xml_tree = etree.ElementTree()
    xml_tree.parse(desc_file)
    descrips = xml_tree.findall("descriptions/desc")

    for descrip in descrips:
        compset_attr = descrip.get("compset")
        if compset_attr in compset:
## start here return the compset_attr and the element value            


###############################################################################
def _main_func(options, work_dir):
###############################################################################

    """Construct compsets html from an XML file."""
        
    # Initialize variables
    all_compsets = dict()
    desc = dict()
    ordered_compsets = collections.OrderedDict()
    html_dict = dict()
    model_version = options.version[0]
    cimeroot = os.environ['CIMEROOT']
    model = get_model()
    srcroot = os.path.dirname(os.path.abspath(get_cime_root()))
    model_config_file = 'config_component_{0}.xml'.format(model)

    # read in all the component config_compsets.xml files
    files = Files()
    components = files.get_components("COMPSETS_SPEC_FILE")   

    # loop through the components and read in the config_compset.xml file
    for component in components:
        # Determine the compsets file for this component
        config_file = files.get_value("COMPSETS_SPEC_FILE", {"component":component})

        expect((config_file),
               "Cannot find any config_compsets.xml file for %s" %component)

        # Check that file exists on disk
        if (os.path.isfile(config_file)):
            compsets = Compsets(config_file)
            help_text, all_compsets = compsets.return_all_values()

            # get the all_compsets sorted
            ordered_compsets = collections.OrderedDict(sorted(all_compsets.items(), key=lambda t: t[0]))

            # load up the html_dict
            html_dict[component] = { 'compsets' : ordered_compsets }

    # load the desc_files list
    desc_files = dict()
    for path, comps in _cime_comps.iteritems():
        for comp in comps:
            desc_files[comp] = os.path.join(cimeroot, 'src/components', path, comp)

    # add to the desc_files list
    for component in components:
        if 'allactive' in component or 'drv' in component:
            desc_files[component] = os.path.join(cimeroot, 'src/drivers/mct/cime_config', model_config_file)
        else:
            desc_files[component] = os.path.join(srcroot, 'components', component, 'cime_config/config_component.xml')

    # get the mode descriptions for every compset and every corresponding config_component.xml file
    for alias, compset in ordered_compsets.iteritems():
        description = ''
        parts = compset.split('_')
        for comp, desc_file in desc_files.iteritems():
            if comp.upper() in parts:
                expect((desc_file),
                       "Cannot find config_component.xml file for %s" %component)
                description += get_desc(compset, desc_file)




                compobj = Component(desc_file)
                description += compobj.get_description(compset).replace(':','<br/>')
                description = description.replace('-----------------------------WARNING ------------------------------------------------',
                                                  '<b>- WARNING -<b/><br/>')
                description = description.replace('-------------------------------------------------------------------------------------',
                                                  '<br/>-----------<br/>')
        html_dict[comp]['desc'] = { alias : description }

    # load up jinja template
    templateLoader = jinja2.FileSystemLoader( searchpath='{0}/templates'.format(work_dir) )
    templateEnv = jinja2.Environment( loader=templateLoader )

    # TODO - get the cesm_version for the CIME root
    tmplFile = 'compdef2html.tmpl'
    template = templateEnv.get_template( tmplFile )
    templateVars = { 'html_dict'     : html_dict,
                     'today'         : _now,
                     'model_version' : model_version,
                     'help'          : help_text }
        
    # render the template
    comp_tmpl = template.render( templateVars )

    # write the output file
    with open( options.htmlfile[0], 'w') as html:
        html.write(comp_tmpl)

    return 0

###############################################################################

if __name__ == "__main__":

    options = commandline_options()
    work_dir = os.path.join(CIMEROOT,"scripts","Tools","xml2html")
    try:
        status = _main_func(options, work_dir)
        sys.exit(status)
    except Exception as error:
        print(str(error))
        sys.exit(1)




