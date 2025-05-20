#!/usr/bin/env python3
"""
CLI Documentation Generator

This script generates a unified Markdown documentation file from multiple Python
scripts that use argparse for command-line argument parsing.

The generator:
1. Scans a directory for Python files
2. Imports each file and looks for ArgumentParser objects
3. Extracts documentation from each parser
4. Combines all documentation into a single Markdown file
"""

import os
import sys
import importlib.util
import inspect
import argparse
import re
from typing import List, Dict, Any, Optional, Tuple


def find_python_files(directory: str) -> List[str]:
    """Find all Python files in the given directory.
    
    Args:
        directory: The directory to search in
        
    Returns:
        A list of paths to Python files
    """
    python_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                python_files.append(os.path.join(root, file))
    return python_files


def import_script(script_path: str) -> Optional[Any]:
    """Import a Python script from its file path.
    
    Args:
        script_path: Path to the Python script
        
    Returns:
        The imported module, or None if import failed
    """
    try:
        # Extract module name from file path
        module_name = os.path.basename(script_path).replace('.py', '')
        
        # Import the module
        spec = importlib.util.spec_from_file_location(module_name, script_path)
        if spec is None or spec.loader is None:
            print(f"Warning: Could not load spec for {script_path}")
            return None
            
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    except Exception as e:
        print(f"Warning: Failed to import {script_path}: {e}")
        return None


def find_parsers_in_module(module: Any) -> List[Tuple[str, argparse.ArgumentParser]]:
    """Find all ArgumentParser objects in a module.
    
    Args:
        module: The imported Python module
        
    Returns:
        A list of (name, parser) tuples
    """
    parsers = []
    
    # Check if there's a main parser defined directly in the module
    if hasattr(module, 'parser') and isinstance(module.parser, argparse.ArgumentParser):
        parsers.append((module.__name__, module.parser))
    
    # Look for functions that create and return parsers
    for name, obj in inspect.getmembers(module):
        # Check functions that might create parsers
        if inspect.isfunction(obj) and ('parse' in name.lower() or 'arg' in name.lower()):
            try:
                result = obj()
                if isinstance(result, argparse.ArgumentParser):
                    parsers.append((f"{module.__name__}.{name}", result))
            except Exception:
                # Skip functions that fail or require arguments
                pass
        
        # Check for parser instances in classes
        elif inspect.isclass(obj):
            for method_name, method in inspect.getmembers(obj, inspect.isfunction):
                if 'parse' in method_name.lower() or 'arg' in method_name.lower():
                    try:
                        # Try to create an instance and call the method
                        instance = obj()
                        result = method(instance)
                        if isinstance(result, argparse.ArgumentParser):
                            parsers.append((f"{module.__name__}.{obj.__name__}.{method_name}", result))
                    except Exception:
                        # Skip methods that fail or require arguments
                        pass
    
    return parsers


def extract_usage_from_help(help_text: str) -> str:
    """Extract the usage part from a help text.
    
    Args:
        help_text: The full help text
        
    Returns:
        The usage section
    """
    usage_match = re.search(r'usage: (.*?)(?:\n\n|\Z)', help_text, re.DOTALL)
    if usage_match:
        return usage_match.group(1).strip()
    return "Usage information not available"


def extract_description(parser: argparse.ArgumentParser) -> str:
    """Extract the description from a parser.
    
    Args:
        parser: The argparse ArgumentParser
        
    Returns:
        The description
    """
    return parser.description or "No description available"


def get_argument_details(action) -> Dict[str, Any]:
    """Extract details from an argparse action.
    
    Args:
        action: The argparse action
        
    Returns:
        A dictionary with argument details
    """
    details = {
        'name': action.dest,
        'help': action.help or "No help available",
        'default': action.default if action.default is not argparse.SUPPRESS else None,
        'required': action.required,
        'type': getattr(action.type, '__name__', str(action.type)) if action.type else None,
    }
    
    if action.choices:
        details['choices'] = list(action.choices)
    
    if action.option_strings:
        details['options'] = action.option_strings
    
    return details


def generate_markdown_for_parser(name: str, parser: argparse.ArgumentParser) -> str:
    """Generate Markdown documentation for a parser.
    
    Args:
        name: The name of the script or function containing the parser
        parser: The argparse ArgumentParser
        
    Returns:
        Markdown documentation for the parser
    """
    help_text = parser.format_help()
    usage = extract_usage_from_help(help_text)
    description = extract_description(parser)
    
    lines = [
        f"## {name}",
        "",
        f"{description}",
        "",
        "### Usage",
        "```",
        f"{usage}",
        "```",
        "",
        "### Arguments",
        ""
    ]
    
    # Group arguments by their group (positional, optional, custom groups)
    grouped_actions = {}
    
    # Handle positional arguments
    positional_args = [action for action in parser._actions 
                      if action.option_strings == [] and action.dest != "help"]
    if positional_args:
        grouped_actions["Positional Arguments"] = positional_args
    
    # Handle optional arguments and other groups
    for group in parser._action_groups:
        if group.title not in ["positional arguments"]:
            group_actions = [action for action in group._group_actions 
                           if action.dest != "help"]  # Skip help action
            if group_actions:
                grouped_actions[group.title.title()] = group_actions
    
    # Generate documentation for each group
    for group_name, actions in grouped_actions.items():
        lines.append(f"#### {group_name}")
        lines.append("")
        
        for action in actions:
            arg_details = get_argument_details(action)
            
            # Format the argument name/options
            if 'options' in arg_details:
                arg_name = ", ".join(f"`{o}`" for o in arg_details['options'])
            else:
                arg_name = f"`{arg_details['name']}`"
            
            # Add type if available
            arg_type = f" ({arg_details['type']})" if arg_details['type'] else ""
            
            # Add required flag if true
            required = " [required]" if arg_details['required'] else ""
            
            lines.append(f"* {arg_name}{arg_type}{required}: {arg_details['help']}")
            
            # Add default value if available
            if arg_details['default'] is not None:
                lines.append(f"  * Default: `{arg_details['default']}`")
            
            # Add choices if available
            if 'choices' in arg_details:
                choices_str = ", ".join(f"`{c}`" for c in arg_details['choices'])
                lines.append(f"  * Choices: {choices_str}")
            
            lines.append("")
        
        lines.append("")
    
    return "\n".join(lines)


def generate_unified_docs(directory: str, output_file: str) -> bool:
    """Generate unified documentation for all scripts in a directory.
    
    Args:
        directory: Directory containing Python scripts
        output_file: Path to output Markdown file
        
    Returns:
        True if successful, False otherwise
    """
    # Find Python files
    python_files = find_python_files(directory)
    if not python_files:
        print(f"No Python files found in {directory}")
        return False
    
    # Import modules and find parsers
    all_parsers = []
    script_parsers = []  # To store parsers from actual scripts, not the generator itself
    
    for script_path in python_files:
        # Skip the documentation generator script itself
        if os.path.basename(script_path) == os.path.basename(sys.argv[0]):
            continue
            
        module = import_script(script_path)
        if module:
            relative_path = os.path.relpath(script_path, directory)
            parsers = find_parsers_in_module(module)
            for name, parser in parsers:
                script_parsers.append((relative_path, name, parser))
    
    # Sort parsers alphabetically by name
    script_parsers.sort(key=lambda x: x[1].lower())
    all_parsers = script_parsers
    
    if not all_parsers:
        print("No argparse parsers found in any scripts")
        return False
    
    # Generate markdown
    lines = [
        "# Command Line Tools Documentation",
        "",
        "This document provides documentation for command line tools.",
        "",
        "## Table of Contents",
        ""
    ]
    
    # Generate table of contents (now sorted alphabetically)
    for script_path, name, _ in all_parsers:
        lines.append(f"* [{name}](#{name.replace('.', '-').lower()})")
    
    lines.append("\n---\n")
    
    # Generate documentation for each parser (now sorted alphabetically)
    for script_path, name, parser in all_parsers:
        lines.append(f"<!-- Source: {script_path} -->")
        lines.append(generate_markdown_for_parser(name, parser))
        lines.append("---\n")
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write("\n".join(lines))
    
    print(f"Documentation written to {output_file}")
    return True


def main():
    """Parse command line arguments and generate documentation."""
    parser = argparse.ArgumentParser(
        description="Generate unified Markdown documentation for CLI scripts with argparse"
    )
    parser.add_argument(
        "directory",
        help="Directory containing Python scripts"
    )
    parser.add_argument(
        "-o", "--output",
        default="cli_documentation.md",
        help="Output Markdown file (default: %(default)s)"
    )
    args = parser.parse_args()
    
    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a directory")
        return 1
    
    success = generate_unified_docs(args.directory, args.output)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
