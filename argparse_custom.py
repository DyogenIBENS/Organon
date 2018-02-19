

"""A custom help formatter that only prints metavar and choices once."""

from argparse import *


class HelpFormatter(HelpFormatter):
    """Replace argparse default metavar formatting: only write it once."""
    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    # My change here compared to argparse.HelpFormatter:
                    # Append only the option string
                    parts.append('%s' % option_string)
                # Now append the METAVAR (or choices)
                parts[-1] += ' %s' % args_string

            return ', '.join(parts)


class RawTextHelpFormatter(HelpFormatter, RawTextHelpFormatter):
    pass

class RawDescriptionHelpFormatter(HelpFormatter, RawDescriptionHelpFormatter):
    pass

class ArgumentDefaultsHelpFormatter(HelpFormatter, ArgumentDefaultsHelpFormatter):
    pass

class ArgumentParser(ArgumentParser):
    def __init__(self, prog=None, usage=None, description=None, epilog=None,
                 parents=[], formatter_class=HelpFormatter, prefix_chars='-',
                 fromfile_prefix_chars=None, argument_default=None,
                 conflict_handler='error', add_help=True, allow_abbrev=True):
        super().__init__(prog, usage, description, epilog, parents,
                         formatter_class, prefix_chars, fromfile_prefix_chars,
                         argument_default, conflict_handler, add_help,
                         allow_abbrev)
