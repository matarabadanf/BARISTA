#!/usr/bin/env python3
"""
Python Code Documentation Generator

This script generates comprehensive Markdown documentation for Python scripts,
including both:
1. Command-line arguments (from argparse)
2. Class structures (methods, properties, docstrings)

The generator scans a directory for Python files, analyzes their contents,
and produces a single unified documentation file.
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
        sys.modules[module_name] = module  # Add to sys.modules to handle imports within the module
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
            
            # Check if the class has a parser as a class attribute
            if hasattr(obj, 'parser') and isinstance(obj.parser, argparse.ArgumentParser):
                parsers.append((f"{module.__name__}.{obj.__name__}.parser", obj.parser))
    
    return parsers


def find_classes_in_module(module: Any) -> List[Tuple[str, type]]:
    """Find all classes in a module.
    
    Args:
        module: The imported Python module
        
    Returns:
        A list of (name, class) tuples
    """
    classes = []
    
    for name, obj in inspect.getmembers(module, inspect.isclass):
        # Skip imported classes (only document classes defined in this module)
        if obj.__module__ == module.__name__:
            classes.append((f"{module.__name__}.{name}", obj))
    
    return classes


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
        f"## CLI Arguments: {name}",
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


def format_docstring(docstring: str) -> str:
    """Format a docstring for Markdown display.
    
    Args:
        docstring: The docstring to format
        
    Returns:
        Formatted docstring for Markdown
    """
    if not docstring:
        return "_No documentation provided_"
    
    # Clean up the docstring
    doc = inspect.cleandoc(docstring)
    
    # Look for parameter descriptions in the docstring
    param_section = False
    formatted_lines = []
    param_descriptions = {}
    current_param = None
    
    for line in doc.split('\n'):
        # Look for common parameter section headers
        if re.match(r'^\s*(Args|Parameters|Arguments):\s*$', line, re.IGNORECASE):
            param_section = True
            continue
        
        # Look for other section headers that would end parameters section
        if param_section and re.match(r'^\s*(Returns|Raises|Yields|Examples):\s*$', line, re.IGNORECASE):
            param_section = False
        
        # Extract parameter descriptions
        if param_section:
            # Match parameters in the form "param_name: description" or "param_name (type): description"
            param_match = re.match(r'^\s*([a-zA-Z0-9_]+)(\s*\([^)]*\))?\s*:\s*(.*)$', line)
            if param_match:
                current_param = param_match.group(1)
                description = param_match.group(3).strip()
                param_descriptions[current_param] = description
            elif current_param and line.strip():
                # Continuation of previous parameter description
                param_descriptions[current_param] += " " + line.strip()
    
    return doc, param_descriptions


def get_method_signature(method) -> str:
    """Get a formatted method signature.
    
    Args:
        method: The method to document
        
    Returns:
        Formatted method signature
    """
    try:
        signature = inspect.signature(method)
        return f"{method.__name__}{signature}"
    except (TypeError, ValueError):
        return f"{method.__name__}(...)"


def get_class_property_info(cls, prop_name: str) -> Dict[str, Any]:
    """Get information about a class property.
    
    Args:
        cls: The class
        prop_name: The property name
        
    Returns:
        Dictionary with property information
    """
    info = {'name': prop_name, 'docstring': None, 'type_hint': None}
    
    # Try to get the property object
    prop = getattr(cls, prop_name, None)
    if isinstance(prop, property):
        # Extract docstring
        info['docstring'] = prop.__doc__
        
        # Try to get type hint from annotations
        if hasattr(cls, '__annotations__') and prop_name in cls.__annotations__:
            info['type_hint'] = cls.__annotations__[prop_name]
    
    return info


def generate_markdown_for_class(name: str, cls: type) -> str:
    """Generate Markdown documentation for a class.
    
    Args:
        name: The full name of the class
        cls: The class object
        
    Returns:
        Markdown documentation for the class
    """
    lines = [f"## Class: {name}", ""]
    
    # Class docstring
    if cls.__doc__:
        doc, _ = format_docstring(cls.__doc__)
        lines.append(doc)
    else:
        lines.append("_No class documentation provided_")
    
    lines.append("")
    
    # Check for base classes
    bases = cls.__bases__
    if bases and not all(base is object for base in bases):
        base_names = [base.__name__ for base in bases if base is not object]
        if base_names:
            lines.append(f"**Inherits from:** {', '.join(base_names)}")
            lines.append("")
    
    # Properties
    properties = []
    for name, member in inspect.getmembers(cls):
        if isinstance(member, property):
            properties.append(name)
    
    if properties:
        lines.append("### Properties")
        lines.append("")
        
        for prop_name in sorted(properties):
            prop_info = get_class_property_info(cls, prop_name)
            type_hint = f" ({prop_info['type_hint']})" if prop_info['type_hint'] else ""
            lines.append(f"#### `{prop_name}`{type_hint}")
            lines.append("")
            
            if prop_info['docstring']:
                doc, _ = format_docstring(prop_info['docstring'])
                lines.append(doc)
            else:
                lines.append("_No documentation provided_")
            
            lines.append("")
    
    # Methods
    methods = []
    for name, member in inspect.getmembers(cls, inspect.isfunction):
        # Skip special methods and private methods (starting with _)
        if not name.startswith('_'):
            methods.append((name, member))
    
    if methods:
        lines.append("### Methods")
        lines.append("")
        
        for method_name, method in sorted(methods, key=lambda x: x[0]):
            signature = get_method_signature(method)
            lines.append(f"#### `{signature}`")
            lines.append("")
            
            if method.__doc__:
                doc, param_descriptions = format_docstring(method.__doc__)
                lines.append(doc)
                
                # Add parameter descriptions if available and not already in docstring
                if param_descriptions and "Args:" not in doc and "Parameters:" not in doc:
                    lines.append("")
                    lines.append("**Parameters:**")
                    lines.append("")
                    for param, desc in param_descriptions.items():
                        lines.append(f"* `{param}`: {desc}")
            else:
                lines.append("_No documentation provided_")
            
            lines.append("")
    
    return "\n".join(lines)


def generate_unified_docs(directory: str, output_file: str, include_classes: bool = True) -> bool:
    """Generate unified documentation for all scripts in a directory.
    
    Args:
        directory: Directory containing Python scripts
        output_file: Path to output Markdown file
        include_classes: Whether to include class documentation
        
    Returns:
        True if successful, False otherwise
    """
    # Find Python files
    python_files = find_python_files(directory)
    if not python_files:
        print(f"No Python files found in {directory}")
        return False
    
    # Import modules and find parsers and classes
    all_parsers = []
    all_classes = []
    
    for script_path in python_files:
        # Skip the documentation generator script itself
        if os.path.basename(script_path) == os.path.basename(sys.argv[0]):
            continue
            
        module = import_script(script_path)
        if module:
            relative_path = os.path.relpath(script_path, directory)
            
            parsers = find_parsers_in_module(module)
            for name, parser in parsers:
                all_parsers.append((relative_path, name, parser))
            
            if include_classes:
                classes = find_classes_in_module(module)
                for name, cls in classes:
                    all_classes.append((relative_path, name, cls))
    
    # Sort parsers and classes alphabetically by name
    all_parsers.sort(key=lambda x: x[1].lower())
    all_classes.sort(key=lambda x: x[1].lower())
    
    if not all_parsers and not all_classes:
        print("No parsers or classes found in any scripts")
        return False
    
    # Generate markdown
    lines = [
        "# Python Code Documentation",
        "",
        "This document provides comprehensive documentation for Python scripts.",
        "",
        "## Table of Contents",
        ""
    ]
    
    # CLI Arguments section in TOC
    if all_parsers:
        lines.append("### Command Line Interfaces")
        lines.append("")
        for script_path, name, _ in all_parsers:
            lines.append(f"* [{name}](#{name.replace('.', '-').lower()})")
        lines.append("")
    
    # Classes section in TOC
    if all_classes:
        lines.append("### Classes")
        lines.append("")
        for script_path, name, _ in all_classes:
            lines.append(f"* [{name}](#{name.replace('.', '-').lower()})")
        lines.append("")
    
    lines.append("---")
    lines.append("")
    
    # Generate CLI Arguments documentation
    if all_parsers:
        lines.append("# Command Line Interfaces")
        lines.append("")
        for script_path, name, parser in all_parsers:
            lines.append(f"<!-- Source: {script_path} -->")
            lines.append(generate_markdown_for_parser(name, parser))
            lines.append("---")
            lines.append("")
    
    # Generate Class documentation
    if all_classes:
        lines.append("# Classes")
        lines.append("")
        for script_path, name, cls in all_classes:
            lines.append(f"<!-- Source: {script_path} -->")
            lines.append(generate_markdown_for_class(name, cls))
            lines.append("---")
            lines.append("")
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write("\n".join(lines))
    
    print(f"Documentation written to {output_file}")
    return True


def main():
    """Parse command line arguments and generate documentation."""
    parser = argparse.ArgumentParser(
        description="Generate unified Markdown documentation for Python scripts"
    )
    parser.add_argument(
        "directory",
        help="Directory containing Python scripts"
    )
    parser.add_argument(
        "-o", "--output",
        default="python_documentation.md",
        help="Output Markdown file (default: %(default)s)"
    )
    parser.add_argument(
        "--no-classes",
        action="store_true",
        help="Skip class documentation"
    )
    parser.add_argument(
        "--cli-only",
        action="store_true",
        help="Only include command line interface documentation"
    )
    args = parser.parse_args()
    
    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a directory")
        return 1
    
    include_classes = not (args.no_classes or args.cli_only)
    success = generate_unified_docs(args.directory, args.output, include_classes)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
